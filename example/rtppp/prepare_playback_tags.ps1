param(
    [string]$RoverPath = 'D:\WORK-FRR\DATA\ROVER',
    [string]$B2bPath = 'D:\WORK-FRR\DATA\ppp',
    [string]$OutputDirectory = $PSScriptRoot,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

function Get-Bits {
    param([byte[]]$Data, [long]$BitPosition, [int]$Length)
    [uint64]$value = 0
    for ($i = 0; $i -lt $Length; $i++) {
        $p = $BitPosition + $i
        $bit = ([int]$Data[[int]($p -shr 3)] -shr (7 - ($p -band 7))) -band 1
        $value = ($value -shl 1) -bor [uint64]$bit
    }
    return $value
}

function Read-RtcmFrames {
    param([string]$Path, [ValidateSet('Rover','B2b')] [string]$Kind)
    $data = [IO.File]::ReadAllBytes($Path)
    $frames = [Collections.Generic.List[object]]::new()
    $offset = 0
    while ($offset -le $data.Length - 6) {
        if ($data[$offset] -ne 0xD3) { $offset++; continue }
        $length = (([int]$data[$offset + 1] -band 3) -shl 8) -bor [int]$data[$offset + 2]
        $end = $offset + $length + 6
        if ($length -lt 2 -or $end -gt $data.Length) { $offset++; continue }
        $type = (([int]$data[$offset + 3]) -shl 4) -bor (([int]$data[$offset + 4]) -shr 4)
        $stamp = $null
        if ($Kind -eq 'Rover' -and $type -ge 1071 -and $type -le 1077) {
            $stamp = [double](Get-Bits $data ($offset * 8L + 48) 30) / 1000.0
        }
        elseif ($Kind -eq 'B2b' -and $type -eq 4047 -and $length -ge 71) {
            # The 4047/64 body is stored as little-endian uint32 words.
            [uint32]$word = (([uint32]$data[$offset + 13]) -shl 24) -bor
                             (([uint32]$data[$offset + 12]) -shl 16) -bor
                             (([uint32]$data[$offset + 11]) -shl 8) -bor
                              [uint32]$data[$offset + 10]
            $outerMessageType = [int]$data[$offset + 9]
            $innerMessageType = [int](($word -shr 26) -band 0x3F)
            if ($outerMessageType -ge 1 -and $outerMessageType -le 4 -and
                $innerMessageType -eq $outerMessageType) {
                $stamp = [double](($word -shr 9) -band 0x1FFFF) + 14.0 # BDT -> GPST
            }
        }
        $frames.Add([pscustomobject]@{ End = [uint32]$end; Stamp = $stamp; Type = $type })
        $offset = $end
    }
    if ($frames.Count -eq 0) { throw "No RTCM frames found in $Path" }
    return ,$frames
}

function Expand-And-UnwrapStamps {
    param([object[]]$Frames, [double]$Period)
    $known = @($Frames | Where-Object { $null -ne $_.Stamp })
    if ($known.Count -eq 0) { throw 'No usable embedded timestamps found.' }
    [double]$previous = $known[0].Stamp
    [double]$wrap = 0
    foreach ($frame in $Frames) {
        if ($null -ne $frame.Stamp) {
            [double]$candidate = $frame.Stamp + $wrap
            if ($candidate -lt $previous - $Period / 2) { $wrap += $Period; $candidate += $Period }
            $frame.Stamp = $candidate
            $previous = $candidate
        }
    }
    [double]$first = ($Frames | Where-Object { $null -ne $_.Stamp } | Select-Object -First 1).Stamp
    [double]$current = $first
    foreach ($frame in $Frames) {
        if ($null -ne $frame.Stamp) { $current = $frame.Stamp }
        else { $frame.Stamp = $current }
    }
}

function Write-TimeTag {
    param([object[]]$Frames, [string]$Path, [double]$BaseStamp,
          [DateTime]$GpsReferenceDate)
    if ((Test-Path -LiteralPath $Path) -and -not $Force) {
        throw "Tag file already exists: $Path (use -Force to replace it)"
    }
    $parent = Split-Path -Parent $Path
    [IO.Directory]::CreateDirectory($parent) | Out-Null
    $stream = [IO.File]::Open($Path, [IO.FileMode]::Create, [IO.FileAccess]::Write)
    $writer = [IO.BinaryWriter]::new($stream)
    try {
        $header = [byte[]]::new(64)
        $label = [Text.Encoding]::ASCII.GetBytes('TIMETAG RTKLIB 2.4.3')
        [Array]::Copy($label, $header, $label.Length)
        # tick_f remains zero in the last four header bytes for both streams.
        $writer.Write($header)
        $epoch = [DateTimeOffset]::new([DateTime]::SpecifyKind($GpsReferenceDate, 'Utc'))
        $writer.Write([uint32]$epoch.ToUnixTimeSeconds())
        $writer.Write([double]0.0)
        foreach ($frame in $Frames) {
            $tick = [Math]::Max(0.0, ($frame.Stamp - $BaseStamp) * 1000.0)
            $writer.Write([uint32][Math]::Round($tick))
            $writer.Write([uint32]$frame.End)
        }
    }
    finally { $writer.Dispose(); $stream.Dispose() }
}

$roverFrames = Read-RtcmFrames $RoverPath Rover
$b2bFrames = Read-RtcmFrames $B2bPath B2b
Expand-And-UnwrapStamps $roverFrames 604800.0

$roverFirst = [double]($roverFrames | Where-Object { $_.Type -ge 1071 -and $_.Type -le 1077 } | Select-Object -First 1).Stamp
$gpsDay = [Math]::Floor($roverFirst / 86400.0)
foreach ($frame in $b2bFrames) { $frame.Stamp += $gpsDay * 86400.0 }
Expand-And-UnwrapStamps $b2bFrames 86400.0

$baseStamp = [Math]::Min([double]$roverFrames[0].Stamp, [double]$b2bFrames[0].Stamp)
$modifiedDate = (Get-Item -LiteralPath $RoverPath).LastWriteTime.Date
$gpsDow = [int]([Math]::Floor($roverFirst / 86400.0) % 7)
$referenceDate = $modifiedDate
for ($delta = -3; $delta -le 3; $delta++) {
    $candidate = $modifiedDate.AddDays($delta)
    if ([int]$candidate.DayOfWeek -eq $gpsDow) { $referenceDate = $candidate; break }
}
$referenceDate = $referenceDate.AddSeconds($baseStamp % 86400.0)

$roverTag = Join-Path $OutputDirectory 'ROVER.tag'
$b2bTag = Join-Path $OutputDirectory 'ppp.tag'
Write-TimeTag $roverFrames $roverTag $baseStamp $referenceDate
Write-TimeTag $b2bFrames $b2bTag $baseStamp $referenceDate

Write-Host "ROVER: $($roverFrames.Count) frames -> $roverTag"
Write-Host "PPP-B2b: $($b2bFrames.Count) frames -> $b2bTag"
Write-Host "Replay reference (GPST): $($referenceDate.ToString('yyyy-MM-dd HH:mm:ss'))"
