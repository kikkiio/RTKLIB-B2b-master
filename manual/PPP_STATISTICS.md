# C++ 定位精度与收敛统计工具

源码：tools/ppp_statistics.cpp，独立 C++17 程序，不依赖 RTKLIB、Python 或第三方库。读取 RTKLIB .pos 文件，输出与示例表格对应的定位质量、数量、百分比、RMS、误差分位数、坐标差分速度及收敛时间。所有距离单位为米，时间单位为秒。

## 本算例统计口径

- 对 260909-2 静态 IF1213 的 full/out/solution.pos 统计，默认选取 Q=6（PPP）计算参考坐标和收敛时刻。
- 以最后 180 秒 PPP 解的 ECEF 均值为固定参考真值，区间为 (末历元−180秒, 末历元]。本例为 04:58:35 至 05:01:34，共 180 个 1 Hz 历元。
- 将位置减去参考 ECEF 后，旋转到参考点的固定 ENU 坐标系。水平 RMS = sqrt(mean(E²+N²))，高程 RMS = sqrt(mean(U²))。这是对固定参考位置的 RMS，不是扣除窗口均值的标准差。
- 在每个 PPP 历元 t 检查向前 300 秒窗口 [t,t+300秒)，本例每窗 300 点，含 t+299 秒、不含 t+300 秒。
- 水平 RMS ≤ 0.10 m 且高程 RMS ≤ 0.20 m 为窗口达标。取从此开始所有后续完整窗口均达标的最早窗口起点为收敛时刻；短暂达标后又超限不能算作收敛。浮点数比较仅使用 1e-12 m 的舍入容差。
- 尾部不足 300 秒的候选窗口不参与判定；最后一个完整窗口仍覆盖数据尾部。收敛后总体 RMS 使用从收敛时刻到末历元的全部数据。
- 采样间隔由目标解时间差的中位数估计。间断超过 1.5 倍采样间隔会使包含该间断的窗口无效，可通过 --max-gap 指定秒数。参考区间必须覆盖足够时长且无这种间断，否则报错。
- 收敛时间主要从采集首个观测历元 2026/09/09 03:40:38.090 GPST 计；同时输出从解算文件首历元 03:41:17 计时的结果。首解前的 38.910 秒不在 .pos 的样本数量内。
- 首达标窗口起点是事后统计的收敛时刻，其窗口完成边界晚 300 秒。确认“后续全部窗口持续达标”需要整段数据。
- 该参考真值按用户指定方法构造，结果衡量相对于末尾均值的内部一致性。

## 编译和运行

在仓库根目录运行：

~~~powershell
cmake --build . --config Release --target ppp_statistics
.\bin\Release\ppp_statistics.exe --self-test
.\example\rtppp\kxw_260909_if1213\run_statistics.ps1
~~~

也可直接调用程序：

~~~powershell
.\bin\Release\ppp_statistics.exe --input "example/rtppp/kxw_260909_if1213/full/out/solution.pos" --output "example/rtppp/kxw_260909_if1213/full/statistics/convergence" --label "260909-2" --start "2026/09/09 03:40:38.090" --reference-seconds 180 --window-seconds 300 --horizontal 0.10 --vertical 0.20
~~~

不提供 --start 时，从输入 .pos 文件首历元计时。--quality 可指定其他 RTKLIB 解质量编号（1..7）。不同质量解不会混入 PPP 的参考均值和收敛判定。输入按时间严格递增，重复时间、无效坐标、缺少参考数据会报错。默认从表头识别 ECEF 或十进制度经纬度（椭球高），也可显式指定 --format ecef/llh。支持空白或逗号字段分隔；暂不支持周/周内秒、度分秒、ENU 基线坐标、其他时区转换。--start 必须与输入时间系统一致。

新建独立构建目录时先执行 cmake -S . -B build；构建命令改为 cmake --build build --config Release --target ppp_statistics。可执行文件仍按工程配置输出到仓库 bin/Release。

在 Visual Studio 开发者命令行中可单文件编译：

~~~text
cl /std:c++17 /EHsc /utf-8 /D_CRT_SECURE_NO_WARNINGS tools\ppp_statistics.cpp /Fe:ppp_statistics.exe
~~~

## 输出文件

前缀由 --output 指定，目录不存在时创建。同名前缀会更新已有统计产物，原始 .pos 不会改动。

| 文件 | 内容 |
|---|---|
| convergence.csv | 类似示例的统计表，分“全时段”和“收敛后”，按解质量列行 |
| convergence.md | 可读报告，包含参考真值、计时起点、判据、首达标和前一窗口 RMS、统计表 |
| convergence_windows.csv | 每个完整候选窗口的起止时间、样本数、有效性、水平/高程 RMS、当前是否达标及以后是否持续达标 |
| convergence_epochs.csv | 各历元 E/N/U 误差、水平模长、高程绝对误差、是否收敛后、是否参与参考均值 |

CSV 使用 UTF-8 BOM，可直接用 Excel 打开。无样本或无法计算时标记 N/A，不将缺失数据伪造为零误差。若没有任何持续达标的完整窗口，报告“未收敛”，收敛时间为 N/A，也不输出“收敛后”统计区间。

数量占比的分母为相应区间中输入 .pos 的全部记录数（含 Q=0 的无效解），不是原始观测总数。水平 CEP68/95/99.7 是水平误差模长的经验分位数，高程 P68/95/99.7 是 |U| 的经验分位数，按排序后的 (n−1)p 位置线性插值。分位数不要求误差服从高斯分布。

坐标差分速度 RMS = sqrt(mean(||(XYZ_i−XYZ_(i−1))/Δt||²))，只使用统计区间内时间相邻、质量相同、未跨越间断的历元对。本 .pos 没有速度字段，因此该指标是坐标变化率，并非接收机多普勒测速精度。表中同时给出实际速度样本对数。

## 验证

内置 --self-test 覆盖已知 ENU 方向、LLH/ECEF 互转、水平向量 RMS、保留固定偏差的高程 RMS、速度差分、分位数、末尾参考窗口边界、数据间断、短数据、临时达标后超限、阈值等号、日期跨天及输入异常。已接入 CTest。

本算例结果另外以独立 NumPy 直接逐窗计算复核；复核信息保存在 full/statistics/verification.json。
