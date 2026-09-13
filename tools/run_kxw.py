"""Replay tagged KXW B2b + RTCM3 data and generate C++ convergence statistics."""
import argparse, hashlib, json, math, re, subprocess, time
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
FORMAT='%Y/%m/%d %H:%M:%S.%f'
def run_json(exe,path):
    p=subprocess.run([str(exe),str(path)],capture_output=True,text=True,check=True)
    return json.loads(p.stdout)
def set_value(text,key,value):
    text,n=re.subn(r'^'+re.escape(key)+r'\s*=[^\r\n]*',lambda m:key+' = '+str(value),text,flags=re.M)
    return text if n else text+'\n'+key+' = '+str(value)+'\n'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',required=True,type=Path,help='Folder containing ROVER, ROVER.tag, PPP, PPP.tag')
    ap.add_argument('--output',required=True,type=Path,help='Fresh result directory')
    ap.add_argument('--template',type=Path,default=ROOT/'example/rtppp/conf/kxw_if1213.conf')
    ap.add_argument('--antenna',type=Path,default=ROOT/'example/rtppp/conf/igs20.atx')
    ap.add_argument('--elevation',type=float,default=7)
    ap.add_argument('--robust',type=int,choices=[0],default=0)
    ap.add_argument('--speed',type=float,default=20)
    ap.add_argument('--trace',type=int,default=2)
    ap.add_argument('--mode',choices=['static','kinematic'],default='static')
    a=ap.parse_args();data=a.data.resolve();out=a.output.resolve();antenna=a.antenna.resolve()
    if not 1<=a.speed<=20 or not 0<=a.elevation<90:ap.error('speed must be 1..20 and elevation 0..<90')
    if out.exists() and any(out.iterdir()):raise RuntimeError('Output directory must be fresh: '+str(out))
    inputs=[data/n for n in ('ROVER','ROVER.tag','PPP','PPP.tag')]
    for p in inputs+[antenna,a.template]:
        if not p.is_file():raise FileNotFoundError(p)
    exe=ROOT/'bin/Release/rtppp.exe';stats=ROOT/'bin/Release/ppp_statistics.exe'
    audit=run_json(ROOT/'bin/Release/inspect_rtcm_input.exe',data/'ROVER')
    correction=run_json(ROOT/'bin/Release/inspect_kxw_input.exe',data/'PPP')
    first=datetime.strptime(audit['first_gpst'],FORMAT);last=datetime.strptime(audit['last_gpst'],FORMAT)
    if audit['decode_errors']:raise RuntimeError('Observation decode errors: '+str(audit))
    if not correction['correction_events']:raise RuntimeError('No KXW B2b corrections decoded')
    (out/'out').mkdir(parents=True)
    text=a.template.resolve().read_text(encoding='utf8')
    settings={'prcopt.mode':8 if a.mode=='static' else 7,'prcopt.elmin':a.elevation,
        'filopt.satantp':antenna.as_posix(),
        'strpath[0]':(data/'ROVER').as_posix()+f'::T::x{a.speed:g}',
        'strpath[2]':(data/'PPP').as_posix()+f'::T::x{a.speed:g}',
        'strpath[3]':(out/'out/solution.pos').as_posix(),'filopt.trace':(out/'run.trace').as_posix(),
        'filopt.tempdir':(out/'temp').as_posix(),
        'prcopt.replay_end':(last+timedelta(seconds=1)).strftime('%Y %m %d %H %M %S')}
    for key,val in settings.items():text=set_value(text,key,val)
    (out/'temp').mkdir();conf=out/'run.conf';conf.write_text(text,encoding='utf8')
    # Allow stream startup and drain; verify data completion independently below.
    seconds=math.ceil((last-first).total_seconds()/a.speed)+60
    command=[str(exe),'-s','-nc','-o',str(conf),'-r','1','-t',str(a.trace),'--run-seconds',str(seconds)]
    print(json.dumps({'output':str(out),'observation_audit':audit,'correction_audit':correction,'run_seconds':seconds}),flush=True)
    start=time.monotonic()
    with (out/'console.log').open('w',encoding='utf8') as log:
        proc=subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=seconds+90,
                            creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    rows=[]
    pos=out/'out/solution.pos'
    if pos.exists():
        for line in pos.read_text().splitlines():
            if re.match(r'^\d{4}/\d{2}/\d{2} ',line):rows.append(line.split())
    epochs=[datetime.strptime(' '.join(r[:2]),FORMAT) for r in rows]
    complete=bool(epochs) and abs((epochs[-1]-last).total_seconds())<=.1
    report={'command':command,'exit_code':proc.returncode,'wall_seconds':time.monotonic()-start,
        'full_input_end_reached':complete,'epochs':len(rows),'quality_counts':dict(Counter(r[5] for r in rows)),
        'first':str(epochs[0]) if epochs else None,'last':str(epochs[-1]) if epochs else None,
        'observation_audit':audit,'correction_audit':correction,'mode':a.mode,
        'input_sha256':{p.name:sha(p) for p in inputs},'antenna_sha256':sha(antenna),
        'executable_sha256':sha(exe),'config_sha256':sha(conf),'solution_sha256':sha(pos) if pos.exists() else None,
        'solution_gaps':[{'from':str(x),'to':str(y),'seconds':(y-x).total_seconds()} for x,y in zip(epochs,epochs[1:]) if (y-x).total_seconds()>1.5]}
    (out/'process.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    if proc.returncode or not complete:raise RuntimeError('Replay did not produce final observation epoch; inspect process.json and console.log')
    if len(epochs)!=len(set(epochs)):raise RuntimeError('Duplicate output epochs')
    if a.mode=='static':
        for h in (.10,.15):
            prefix=out/'statistics'/f'h{round(h*100)}'
            p=subprocess.run([str(stats),'--input',str(pos),'--output',str(prefix),'--label',data.name,
                '--start',audit['first_gpst'],'--horizontal',str(h),'--vertical','.20'],capture_output=True,text=True)
            if p.returncode:raise RuntimeError(p.stderr)
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
if __name__=='__main__':main()
