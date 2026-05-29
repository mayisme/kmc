#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path


def main() -> None:
    out_dir = Path('/home/xiao/文档/KMC/trilayer_barrier_analysis')
    out_dir.mkdir(parents=True, exist_ok=True)

    # Values read from the user-provided schematic.
    wvtr_single_mg_m2_day = 1.26
    wvtr_trilayer_mg_m2_day = 0.046
    d_al2o3 = 5.26e-19
    d_polymer = 2.0e-13
    l_al = 44e-9
    l_poly = 208e-9
    l_total = 2.0 * l_al + l_poly
    defect_diameter = 1e-6
    defect_radius = defect_diameter / 2.0
    defect_spacing = 100e-6

    improvement = wvtr_single_mg_m2_day / wvtr_trilayer_mg_m2_day
    reduction_percent = (1.0 - wvtr_trilayer_mg_m2_day / wvtr_single_mg_m2_day) * 100.0

    # 1D resistance model: R = L / D. This ignores defect misalignment.
    r_single = l_al / d_al2o3
    r_trilayer_series = 2.0 * l_al / d_al2o3 + l_poly / d_polymer
    wvtr_series_mg_m2_day = wvtr_single_mg_m2_day * (r_single / r_trilayer_series)
    series_improvement = wvtr_single_mg_m2_day / wvtr_series_mg_m2_day

    # Equivalent effective diffusivity inferred from WVTR ratio, assuming same boundary concentration
    # and using J = D_eff * DeltaC / L.
    d_eff_trilayer_from_wvtr = (
        d_al2o3
        * (wvtr_trilayer_mg_m2_day / wvtr_single_mg_m2_day)
        * (l_total / l_al)
    )
    d_eff_series = l_total / r_trilayer_series
    extra_factor_beyond_series = wvtr_series_mg_m2_day / wvtr_trilayer_mg_m2_day

    # Defect geometry.
    defect_area_fraction = math.pi * defect_radius**2 / defect_spacing**2
    direct_alignment_probability = defect_area_fraction
    inverse_alignment_probability = 1.0 / direct_alignment_probability

    # Characteristic diffusion times.
    t_al_cross = l_al**2 / (2.0 * d_al2o3)
    t_two_al_cross = 2.0 * t_al_cross
    mean_distance_to_nearest_defect = defect_spacing * (
        math.sqrt(2.0) + math.log(1.0 + math.sqrt(2.0))
    ) / 6.0
    rms_distance_to_nearest_defect = defect_spacing / math.sqrt(6.0)
    t_lateral_mean_distance = mean_distance_to_nearest_defect**2 / (4.0 * d_polymer)
    t_lateral_rms_distance = rms_distance_to_nearest_defect**2 / (4.0 * d_polymer)

    # 2D narrow-escape estimate: time to find a small absorbing target in one defect-spacing cell.
    # Order-one constants depend on boundary conditions, so the result should be used as a scale.
    cell_area = defect_spacing**2
    cell_radius_equiv = defect_spacing / math.sqrt(math.pi)
    narrow_escape_2pi = cell_area / (2.0 * math.pi * d_polymer) * math.log(cell_radius_equiv / defect_radius)
    narrow_escape_4pi = cell_area / (4.0 * math.pi * d_polymer) * math.log(cell_radius_equiv / defect_radius)

    # If WVTR is inversely proportional to a characteristic passage time, infer the search penalty.
    inferred_trilayer_time = improvement * t_al_cross
    inferred_extra_search_time = max(0.0, inferred_trilayer_time - t_two_al_cross)

    results = {
        'input': {
            'WVTR_single_44nm_Al2O3_mg_m2_day': wvtr_single_mg_m2_day,
            'WVTR_trilayer_mg_m2_day': wvtr_trilayer_mg_m2_day,
            'D_Al2O3_m2_s': d_al2o3,
            'D_polymer_m2_s': d_polymer,
            'Al2O3_layer_thickness_m': l_al,
            'polymer_layer_thickness_m': l_poly,
            'defect_diameter_m': defect_diameter,
            'defect_spacing_m': defect_spacing,
        },
        'barrier_effect': {
            'WVTR_reduction_factor_single_over_trilayer': improvement,
            'WVTR_reduction_percent': reduction_percent,
            'trilayer_over_single_WVTR_ratio': wvtr_trilayer_mg_m2_day / wvtr_single_mg_m2_day,
        },
        'one_dimensional_series_model': {
            'single_resistance_L_over_D': r_single,
            'trilayer_resistance_L_over_D': r_trilayer_series,
            'predicted_trilayer_WVTR_mg_m2_day': wvtr_series_mg_m2_day,
            'predicted_reduction_factor': series_improvement,
            'extra_reduction_factor_beyond_series': extra_factor_beyond_series,
            'D_eff_series_m2_s': d_eff_series,
            'D_eff_inferred_from_observed_WVTR_m2_s': d_eff_trilayer_from_wvtr,
        },
        'defect_misalignment_model': {
            'defect_area_fraction': defect_area_fraction,
            'direct_vertical_alignment_probability': direct_alignment_probability,
            'inverse_alignment_probability': inverse_alignment_probability,
            'mean_distance_to_nearest_defect_m': mean_distance_to_nearest_defect,
            'rms_distance_to_nearest_defect_m': rms_distance_to_nearest_defect,
            'lateral_search_time_mean_distance_s': t_lateral_mean_distance,
            'lateral_search_time_rms_distance_s': t_lateral_rms_distance,
            'narrow_escape_search_time_2pi_s': narrow_escape_2pi,
            'narrow_escape_search_time_4pi_s': narrow_escape_4pi,
            'Al2O3_single_crossing_time_s': t_al_cross,
            'two_Al2O3_crossing_time_s': t_two_al_cross,
            'inferred_trilayer_characteristic_time_s': inferred_trilayer_time,
            'inferred_extra_search_time_s': inferred_extra_search_time,
        },
    }

    (out_dir / 'trilayer_barrier_analysis.json').write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8'
    )

    md = f"""# 三层膜阻滞效果计算

## 输入

- 单层 44 nm Al2O3 的 WVTR: {wvtr_single_mg_m2_day:g} mg m^-2 day^-1
- 三层膜 WVTR: {wvtr_trilayer_mg_m2_day:g} mg m^-2 day^-1
- D_Al2O3 = {d_al2o3:.3e} m^2/s
- D_polymer = {d_polymer:.3e} m^2/s
- Al2O3 / polymer / Al2O3 = 44 nm / 208 nm / 44 nm
- Al2O3 缺陷直径 1 um，缺陷间距 100 um

## 直接阻滞效果

- WVTR 降低倍数 = {improvement:.3f}x
- WVTR 降低百分比 = {reduction_percent:.2f}%

## 一维串联扩散只能解释一小部分

用 R = L / D 估计串联阻力：

- 单层 Al2O3 阻力 R_single = {r_single:.3e}
- 三层串联阻力 R_trilayer = {r_trilayer_series:.3e}
- 预测三层 WVTR = {wvtr_series_mg_m2_day:.3f} mg m^-2 day^-1
- 串联扩散只能给出约 {series_improvement:.2f}x 降低

实验/图中三层膜实际是 {improvement:.2f}x 降低，因此还需要一个额外的
{extra_factor_beyond_series:.2f}x 阻滞因子。这个因子不能由厚度串联解释，
应主要来自缺陷错位后的横向搜索、窄逃逸和路径曲折。

## 缺陷错位与二维搜索

- 缺陷面积分数 f = pi(0.5 um)^2 / (100 um)^2 = {defect_area_fraction:.3e}
- 若上下缺陷完全随机，直接垂直对齐概率约为 {direct_alignment_probability:.3e}
- 1/f = {inverse_alignment_probability:.0f}，说明直接贯通孔非常罕见

有机层 D 很大，但水分子不是只跨过 208 nm 厚度，而是要在二维平面中寻找
下一个 1 um 尺度缺陷入口。特征横向搜索时间：

- 平均最近缺陷距离约 {mean_distance_to_nearest_defect * 1e6:.1f} um，t ~ r^2/(4D) = {t_lateral_mean_distance:.0f} s
- RMS 距离约 {rms_distance_to_nearest_defect * 1e6:.1f} um，t ~ {t_lateral_rms_distance:.0f} s
- 2D 窄逃逸估计 t ~ A/(2piD) ln(R/a) = {narrow_escape_2pi:.0f} s
- 另一常用常数口径 A/(4piD) ln(R/a) = {narrow_escape_4pi:.0f} s

由 WVTR 反推的额外搜索时间约 {inferred_extra_search_time:.0f} s，
与二维窄逃逸估计同量级，因此“缺陷错位 + 横向搜索”可以解释三层膜远强于
简单串联扩散的阻滞效果。

## 建议的计算分析路线

1. 先用 WVTR 比值给出宏观阻滞倍数。
2. 用 R = L/D 做一维串联扩散基线，证明厚度和 D 的串联只能解释约 2x。
3. 引入缺陷几何：缺陷直径、间距、面积分数和上下缺陷对齐概率。
4. 用二维扩散/窄逃逸模型计算有机层中的横向搜索时间。
5. 若要更贴近真实结构，用 kMC 或 Brownian dynamics 生成路径分布、FPT、Ptrans、Lpath、tau、D_eff 和 WVTR。
"""
    (out_dir / 'trilayer_barrier_analysis.md').write_text(md, encoding='utf-8')
    print(json.dumps(results, ensure_ascii=False, indent=2))
    print(f'outputs: {out_dir}')


if __name__ == '__main__':
    main()
