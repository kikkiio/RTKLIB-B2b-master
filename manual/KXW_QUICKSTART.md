> 此快照尚未加入 Huber 功能，只支持 --robust 0；下列残差开关的新增功能说明适用于第三版。

# KXW 数据自行解算

输入目录需要 ROVER、ROVER.tag、PPP、PPP.tag，分别是 RTCM3 观测及 KXW B2b 改正与配套时间标签。不要修改原始文件名，也不要只复制主文件漏掉 .tag。卫星 ANTEX 需覆盖采集日期和所用频段；接收机型号未知的功能保留。

## Windows 构建

在仓库根目录执行，需要 Visual Studio C++ 工具链、CMake、Python 3（调用脚本只使用标准库）：

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target rtppp inspect_rtcm_input inspect_kxw_input ppp_statistics --parallel 4
```

程序输出在仓库根目录的 bin/Release。本机已有构建可直接运行。

## 一条命令解算并统计

```powershell
python tools/run_kxw.py --data D:/WORK-FRR/DATA/260909-4 --output example/rtppp/my_run_260909_4 --antenna example/rtppp/conf/igs20_2350.atx --elevation 3 --robust 0
```

- 换数据：只修改 `--data`。
- 换结果目录：修改 `--output`，必须是新的空目录，防止旧滤波状态和输出混入。
- 换高度角：`--elevation 3` 或 `--elevation 10`。
- 原有剔除：`--robust 0`；新增降权加原有剔除：`--robust 0`。
- 默认静态：`--mode static`。确为动态数据时用 `--mode kinematic`；动态模式不自动生成“末尾均值为真值”的静态统计。
- 默认20倍回放：处理能力不足时用 `--speed 1`、`5` 或 `10`，支持1至20。
- 其他模型参数：编辑 `example/rtppp/conf/kxw_if1213.conf`，或用 `--template 自己的.conf`。脚本仅覆盖数据/结果路径、卫星天线文件、模式、高度角、降权参数和回放终点，保留模板中的频率组合、接收机天线设置等。

脚本先用生产解码器读取日期及首末观测历元，自动设置回放终点和运行时长；回放结束后核对是否输出到最后观测历元。异常结束或未到达末历元会报错，不能只凭程序退出就判断完成。

## 输出

- `run.conf`：这一次完整配置，可检查实际生效设置。
- `out/solution.pos`：逐历元定位结果，ECEF米，GPST时间。
- `process.json`：输入审计、文件哈希、首末历元、质量数量及输出间断。
- `statistics/h10.md`、`.csv`：水平10 cm/高程20 cm收敛与精度。
- `statistics/h15.md`、`.csv`：水平15 cm/高程20 cm收敛与精度。
- `statistics/*_epochs.csv`：逐历元 ENU；`*_windows.csv`：逐窗口 RMS 与判定。
- `*.trace`、`*.stat`：残差控制、参与卫星和滤波诊断。

统计采用最后180秒目标PPP解的ECEF均值作固定参考，前向300秒RMS窗口持续达标；跨越大于1.5倍中位采样间隔的窗口无效。计时从原始首观测开始。若未收敛，表格明确标注；不会补造结果。自身尾段均值只适合内部一致性评估。

## 直接调用 C/C++ 程序

对脚本已生成的 run.conf，可以自行改路径、时间后直接运行：

```powershell
bin/Release/rtppp.exe -s -nc -o 完整路径/run.conf -r 1 -t 2 --run-seconds 183
bin/Release/ppp_statistics.exe --input 完整路径/out/solution.pos --output 完整路径/statistics/h10 --start "2026/09/09 08:50:19.000" --horizontal 0.10 --vertical 0.20
```

183秒仅适用于本次约2443秒数据、20倍回放加60秒启动/排空余量。直接运行时应从新的结果目录启动，配置中 ROVER/PPP、天线、输出和 trace 使用完整路径，并设置与新数据相符的 `prcopt.replay_end`；否则旧日期可能提前截断新数据。优先使用通用脚本自动处理这些参数。

接收机型号未知保持 `prcopt.anttype[0]` 为空、`prcopt.posopt[1]=0`；已知型号时在模板中设置 ANTEX、型号和几何偏心。GPS/BDS IF 配对分别由 `prcopt.gps_if_pairs`、`prcopt.bds_if_pairs` 配置。
