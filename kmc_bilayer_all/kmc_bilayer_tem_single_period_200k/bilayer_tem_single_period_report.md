# Al2O3/ZnO 双层单周期 kMC 计算报告

## 结构定义

本报告只讨论一个双层周期：

```text
湿侧
|
Al2O3 4.5 nm: 致密非晶 + 稀疏针孔/非贯通自由体积缺陷
|
Al2O3/ZnO interface: 横向搜索层
|
ZnO 6.0 nm: 单层截断纤锌矿晶粒，晶界贯穿 ZnO 厚度
|
干侧
```

该结构对应 TEM 中 5 nm 量级 Al2O3/ZnO 纳米叠层的物理图像：ZnO 层厚度小于或接近单个 ZnO 纳米晶粒尺度，因此 ZnO 不按厚膜式多层晶粒堆叠计算，而按被上下界面截断的单层多晶纤锌矿层计算。ZnO 晶粒本体不可渗透，水汽主要沿晶界、界面错配区域和局部缺陷迁移。

## 传输路径

水分子的计算路径定义为：

```text
表面扩散找 Al2O3 针孔
→ 穿过 Al2O3 针孔
→ 在 Al2O3/ZnO 界面横向搜索 ZnO 晶界入口
→ 沿 ZnO 晶界穿过 6 nm ZnO
→ 到达干侧
```

若粒子没有找到 Al2O3 针孔、没有从界面进入 ZnO 晶界，或进入被阻断的晶界段，则在有限事件数或有限时间内记为未穿透。

## 几何与物理假设

- Al2O3 厚度 = 4.5 nm
- ZnO 厚度 = 6.0 nm
- 总厚度 = 10.5 nm
- 网格分辨率 = 0.5 nm
- 横向宽度 = 320 nm，周期边界
- 粒子数 = 200000
- 最大事件数 = 5000000 / particle
- 温度 = 311.15 K
- 相对湿度 = 0.90
- Al2O3 基体 = 不可渗透
- ZnO 晶粒本体 = 不可渗透
- Al2O3 针孔 = 1-2 nm，间距 50-80 nm
- Al2O3 非贯通自由体积缺陷 = 局部缺陷，不计为表面贯通入口
- ZnO 晶粒模型 = 单层截断横向纤锌矿晶粒
- ZnO 晶界 = 贯穿 ZnO 厚度，但带轻微倾斜/弯曲
- ZnO 晶界阻断比例 = 0.20

## 迁移能垒

- Al2O3 表面扩散 = 0.45 eV
- Al2O3 针孔迁移 = 0.58 eV
- Al2O3/ZnO 界面迁移 = 0.62 eV
- ZnO 晶界迁移 = 0.72 eV
- Al2O3 致密基体 = 不可渗透
- ZnO 晶粒本体 = 不可渗透
- 阻断晶界 = 不可渗透

## 核心计算量

本次输出以下指标：

- Ptrans(t)：给定模拟窗口内的穿透概率。
- FPT：首次通过时间分布。
- Lpath：水分子路径长度。
- tau = Lpath / Lfilm：曲折度。
- D_eff = L_total^2 / (2 * mean_FPT)：有效扩散系数。
- WVTR = D_eff * delta_C / L_total：水汽透过率。
- 未穿透粒子的最终位置分布：用于判断粒子主要滞留在 Al2O3、界面还是 ZnO。

## 计算结果

- Ptrans = 0.023115
- 穿透粒子数 = 4623 / 200000
- FPT mean = 4526.63 s
- FPT median = 4485.18 s
- FPT p10 / p90 = 950.76 / 8175.06 s
- Lpath mean = 1.2465 mm
- Lpath median = 1.2316 mm
- tau mean = 1.187e5
- tau median = 1.173e5
- D_eff, mean-FPT = 1.218e-20 m2/s
- D_eff, median-FPT = 1.229e-20 m2/s
- WVTR, ideal vapor delta_c = 4.164e-06 g m^-2 day^-1
- WVTR, literature sorbed C1 = 6.784e-04 g m^-2 day^-1
- 未穿透粒子数 = 195377
- 未穿透粒子最终深度 mean = 2.570 nm
- 未穿透粒子最终深度 median = 3.000 nm

## 结果解释

本次模型中，ZnO 层不是多层晶粒堆叠，而是一个 6 nm 厚的单层截断纤锌矿晶粒层。晶界可以贯穿 ZnO 厚度，但其入口位置与 Al2O3 针孔出口通常不对齐，因此水分子必须在 Al2O3/ZnO 界面进行横向搜索。该界面搜索与 ZnO 晶界阻断共同决定穿透概率。

Ptrans 仅为 0.023115，说明大部分粒子在有限事件数内没有完成穿透。未穿透粒子的最终深度中位数约为 3.0 nm，表明主要滞留发生在 Al2O3 针孔内部或 Al2O3/ZnO 界面附近。能够穿透的粒子平均路径长度达到 1.2465 mm，相对于 10.5 nm 膜厚的曲折度约为 1.187e5，说明成功路径并非近似直通，而是由表面搜索、界面横向搜索和 ZnO 晶界受限扩散共同拉长。

因此，对 4.5 nm Al2O3 / 6 nm ZnO 双层单周期结构，合理的计算结论是：ZnO 超薄层不能提供厚膜式多层晶粒迷宫，但 Al2O3 针孔与 ZnO 晶界入口的错配、ZnO 晶粒本体不可渗透、以及局部晶界阻断，仍能显著降低有限时间穿透概率并保持较高曲折度。

## 输出文件

- `bilayer_corrected_results.json`: 主 kMC 结果
- `bilayer_corrected_results.csv`: 主 kMC 结果表
- `bilayer_corrected_raw.npz`: 原始粒子数组和几何网格
- `bilayer_corrected_summary.png`: 结构与 FPT 摘要图
- `transport_metrics/kmc_transport_metrics.json`: 后处理指标
- `transport_metrics/kmc_transport_metrics.csv`: 后处理指标表
- `transport_metrics/kmc_metric_distributions.png`: FPT、路径、曲折度与未穿透深度分布
- `transport_metrics/kmc_fpt_vs_path_length.png`: FPT 与路径长度散点图

## 复现命令

```bash
python kmc_bilayer_all/kmc_bilayer_corrected.py \
  --particles 200000 \
  --max-events 5000000 \
  --width-nm 320 \
  --periods 1 \
  --threads-per-block 256 \
  --out-dir kmc_bilayer_all/kmc_bilayer_tem_single_period_200k \
  --seed 20260529

python kmc_metrics_postprocess.py \
  kmc_bilayer_all/kmc_bilayer_tem_single_period_200k/bilayer_corrected_raw.npz \
  --out-dir kmc_bilayer_all/kmc_bilayer_tem_single_period_200k/transport_metrics \
  --film-thickness-nm 10.5 \
  --dx-nm 0.5
```
