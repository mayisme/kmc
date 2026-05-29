# Al2O3/ZnO 双层膜单层截断 ZnO 晶粒 kMC 模型

## 文献约束后的建模思想

ZnO 常见稳定晶体结构为六方纤锌矿。对 ALD 或类似气相沉积 ZnO 薄膜，文献中常见描述是多晶纤锌矿薄膜，晶粒可表现为平行基底的楔形晶粒和垂直基底的细柱状晶粒共存。由于本模型中的 ZnO 层厚度只有 6 nm，而典型 ZnO 晶粒横向尺度约为 10 nm 量级，因此不宜把 ZnO 层画成或算成多个上下堆叠晶粒层。

本次重新计算采用更保守的几何假设：6 nm ZnO 是被膜厚截断的单层横向多晶纤锌矿晶粒域。每个 ZnO 晶粒在厚度方向近似贯穿 6 nm，晶界主要是横向相邻晶粒之间的边界。水分子不能穿过致密 Al2O3 基体和 ZnO 晶粒本体，只能通过 Al2O3 针孔、Al2O3/ZnO 界面层和 ZnO 晶界网络迁移。

## 本次模型

- 结构: Al2O3 / ZnO = 4.5 / 6.0 nm
- 总厚度 = 10.5 nm
- 网格分辨率 = 0.5 nm
- 横向模拟宽度 = 320 nm，周期边界
- ZnO 晶体结构假设 = polycrystalline hexagonal wurtzite ZnO
- ZnO 晶粒模型 = single truncated lateral wurtzite grain layer
- ZnO 平均晶粒尺寸 = 7.5 nm
- ZnO 厚度 / 平均晶粒尺寸 = 0.8
- Al2O3 针孔直径 = 1-2 nm
- Al2O3 针孔间距 = 50-80 nm
- ZnO 晶界阻断比例 = 0.20
- 温度 = 311.15 K，RH = 0.90

## 迁移能垒

- Al2O3 表面扩散: 0.45 eV
- Al2O3 针孔: 0.58 eV
- Al2O3/ZnO 界面层: 0.62 eV
- ZnO 晶界: 0.72 eV
- Al2O3 致密基体: 不可渗透
- ZnO 晶粒本体: 不可渗透
- 阻断晶界: 不可渗透

## 输出

- 粒子数 = 200000
- 穿透粒子数 = 93717
- Ptrans = 0.468585
- FPT mean = 4097 s
- FPT median = 3815 s
- FPT p10 / p90 = 692 / 7997 s
- Lpath mean = 1.100 mm
- Lpath median = 1.025 mm
- tau mean = 1.048e5
- tau median = 9.767e4
- D_eff = 1.345e-20 m2/s
- WVTR, ideal vapor delta_c = 4.600e-06 g m^-2 day^-1
- WVTR, literature sorbed C1 = 7.495e-04 g m^-2 day^-1

## 解释

与旧的二维 Voronoi ZnO 网络相比，新的 ZnO 层不再表示成多个上下堆叠晶粒，而是表示为一个 6 nm 厚的横向多晶层。这个处理更符合“ZnO 膜厚小于或接近单个晶粒尺度”的情况：晶粒被上下表面截断，晶界在厚度方向可贯穿薄层，但横向位置与 Al2O3 针孔不必对齐。

双层膜的主要阻滞仍来自缺陷错配和晶界受限扩散，而不是 ZnO 本体扩散。水分子需要先在 Al2O3 表面找到针孔，穿过 Al2O3 后在界面层横向寻找 ZnO 晶界入口，再沿单层截断 ZnO 晶界向干侧迁移。ZnO 晶粒本体不可渗透，使传输被限制到少数晶界通道；20% 晶界阻断项代表死端、闭合或被局部致密化堵塞的晶界段。

本次 Ptrans 高于旧版 0.1833，原因不是阻隔机理变弱的简单结论，而是几何模型发生了变化：单层截断晶粒模型中，部分 ZnO 晶界在厚度方向贯穿薄层，成功进入晶界的粒子更容易完成穿透；但成功路径仍然有约 1.1 mm 的平均晶格路径长度，曲折度仍达到 1e5 量级。因此，该模型更适合描述 6 nm ZnO 薄层的“薄层晶界贯通 + 横向入口错配”机制，而不是多层晶粒迷宫机制。

## 结果判读

这次重算最重要的变化是几何解释，而不是能垒体系。旧版二维 Voronoi 模型默认 ZnO 内部存在较多上下错开的晶界交叉和死端，适合较厚、多晶粒层叠更明显的 ZnO；新版模型把 6 nm ZnO 视为单层截断晶粒，晶界可在膜厚方向贯通，因此有限时间内成功穿透的粒子比例从 0.1833 升至 0.468585。

但是，成功粒子的平均 FPT 从旧版约 4258 s 变为 4097 s，D_eff 从约 1.29e-20 m2/s 变为 1.35e-20 m2/s，仍处于同一量级。这说明在当前能垒和步长设定下，真正控制 D_eff 的不是“有没有多层 ZnO 晶粒堆叠”，而是 Al2O3 针孔出口、界面层横向搜索和 ZnO 晶界低维通道中的长路径扩散。换言之，6 nm ZnO 不能提供厚膜式多层晶界迷宫，但仍可以通过晶界入口错配和晶粒本体不可渗透来维持较高曲折度。

因此，推荐在后续论文或汇报中采用以下表述：对于 4.5 nm Al2O3 / 6 nm ZnO 双层结构，ZnO 更合理地表示为单层截断的多晶纤锌矿层；其阻隔贡献不是来自多个 ZnO 晶粒层逐层串联，而是来自 Al2O3 针孔与 ZnO 晶界入口的不对齐、界面横向搜索，以及 ZnO 晶粒本体对水汽的排斥。

## 与三层模型的关系

三层 Al2O3/polymer/Al2O3 模型把 polymer 层看作横向 Brownian 搜索惩罚，并用目标 WVTR 校准搜索时间。这里的 Al2O3/ZnO 双层仍只计算两层，但借鉴了三层模型的分析口径：重点不是简单相加厚度，而是追踪分子从上层缺陷到下层可通道之间的横向搜索、路径增长、FPT 增大和 D_eff 降低。

差异在于，本次双层计算没有把 ZnO 层替换成经验校准的 Brownian 搜索时间，而是保留显式网格 kMC。ZnO 晶粒和晶界由几何生成，粒子按 Arrhenius 跳跃速率在允许通道内运动。

## 局限性

- ZnO 晶粒尺度采用 7.5 nm 均值，是对 5-10 nm 量级晶粒的模型化取值；真实晶粒尺寸应由 TEM、AFM 或 XRD Scherrer 分析约束。
- 6 nm ZnO 是否已经充分结晶取决于沉积温度、前驱体、基底和退火条件；若实验显示超薄层部分非晶，应引入非晶 ZnO 或混合相通道。
- 当前 D_eff 主要基于成功穿透粒子的平均 FPT；对未穿透粒子仍属于删失数据，后续应加入有限时间通量或生存分析口径。
- 能垒尚未用目标 WVTR 反演校准，适合作为机制比较和参数敏感性分析，不宜直接当作唯一实验预测值。

## 参考依据

- Britannica 的 zinc oxide 条目指出 ZnO 结晶于 wurtzite 结构，键合兼具离子性和共价性: https://www.britannica.com/science/zinc-oxide
- ZnO 晶体结构综述资料给出六方纤锌矿 ZnO 和晶格参数量级: https://pmc.ncbi.nlm.nih.gov/articles/PMC5109597/
- RSC Journal of Materials Chemistry C 的 ALD ZnO 生长研究显示衍射图谱匹配六方纤锌矿 ZnO，并观察到平行基底的楔形晶粒和垂直基底的细柱状晶粒共存: https://pubs.rsc.org/en/content/articlehtml/2021/tc/d0tc05439a
- ZnO/AlN 纤锌矿薄膜微结构研究显示 ZnO 薄膜可呈 [0001] 织构和柱状微结构，说明薄膜晶界形貌强烈依赖沉积过程和厚度: https://www.sciencedirect.com/science/article/abs/pii/S0955221998004531

## 改进方向

1. 用实验或文献的 ZnO 晶粒尺寸分布替代固定 5-10 nm 范围。
2. 增加多 seed 几何重复，报告 Ptrans、FPT、D_eff 和 WVTR 的均值与置信区间。
3. 对 ZnO 晶界阻断率、界面能垒、晶界能垒、针孔密度做参数扫描。
4. 增加删失校正：同时报告成功粒子 FPT、有限时间 Ptrans(t) 和基于通量的 WVTR。
5. 若有实测 WVTR，用三层模型类似的方式反推有效界面搜索惩罚或校准能垒。

## 输出文件

- `bilayer_corrected_results.json`: 主 kMC 输出
- `bilayer_corrected_raw.npz`: 原始粒子 FPT、步数、穿透标记、最终深度和几何网格
- `bilayer_corrected_summary.png`: 结构与 FPT 摘要图
- `transport_metrics/kmc_transport_metrics.json`: 后处理传输指标
- `transport_metrics/kmc_metric_distributions.png`: FPT、路径长度、曲折度和未穿透深度分布
- `transport_metrics/kmc_fpt_vs_path_length.png`: FPT 与路径长度关系图

## 复现命令

```bash
python kmc_bilayer_all/kmc_bilayer_corrected.py \
  --particles 200000 \
  --max-events 5000000 \
  --width-nm 320 \
  --periods 1 \
  --threads-per-block 256 \
  --out-dir kmc_bilayer_all/kmc_bilayer_single_layer_200k \
  --seed 20260529

python kmc_metrics_postprocess.py \
  kmc_bilayer_all/kmc_bilayer_single_layer_200k/bilayer_corrected_raw.npz \
  --out-dir kmc_bilayer_all/kmc_bilayer_single_layer_200k/transport_metrics \
  --film-thickness-nm 10.5 \
  --dx-nm 0.5
```
