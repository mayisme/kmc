#!/usr/bin/env python3
"""Coarse-grained Brownian/kMC model for Al2O3/polymer/Al2O3 trilayer barriers.

The model is intentionally mesoscopic:

1. Dense Al2O3 is treated as impermeable except at micron-scale defects.
2. A molecule crossing an Al2O3 defect has a 1D first-passage time scale
   L^2/(2D_Al2O3).
3. The polymer interlayer decouples top and bottom Al2O3 defects. Transport
   through the polymer is represented as a lateral Brownian search for the next
   lower Al2O3 defect, calibrated to the observed trilayer WVTR.
4. Path length is reported at a coarse 1 um Brownian-hop resolution. Like any
   Brownian arc length, it is resolution-dependent; FPT, D_eff and WVTR are the
   primary physically stable outputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def saturation_vapor_pressure_pa(temp_c: float) -> float:
    return 611.21 * math.exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))


def water_vapor_density_kg_m3(temp_c: float, rh: float) -> float:
    p = saturation_vapor_pressure_pa(temp_c) * rh
    return p * 0.01801528 / (8.314462618 * (temp_c + 273.15))


def finite_summary(a: np.ndarray) -> dict[str, float | int | None]:
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {
            'count': 0, 'mean': None, 'median': None, 'std': None,
            'p01': None, 'p05': None, 'p10': None, 'p25': None,
            'p75': None, 'p90': None, 'p95': None, 'p99': None,
            'min': None, 'max': None,
        }
    return {
        'count': int(a.size),
        'mean': float(np.mean(a)),
        'median': float(np.median(a)),
        'std': float(np.std(a)),
        'p01': float(np.percentile(a, 1)),
        'p05': float(np.percentile(a, 5)),
        'p10': float(np.percentile(a, 10)),
        'p25': float(np.percentile(a, 25)),
        'p75': float(np.percentile(a, 75)),
        'p90': float(np.percentile(a, 90)),
        'p95': float(np.percentile(a, 95)),
        'p99': float(np.percentile(a, 99)),
        'min': float(np.min(a)),
        'max': float(np.max(a)),
    }


def flatten(prefix: str, obj, rows: list[tuple[str, object]]) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            flatten(f'{prefix}.{k}' if prefix else k, v, rows)
    elif isinstance(obj, list):
        rows.append((prefix, json.dumps(obj, ensure_ascii=False)))
    else:
        rows.append((prefix, obj))


def make_plots(out_dir: Path, arrays: dict[str, np.ndarray], observation_time_s: float) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fpt = arrays['first_passage_time_s']
    lpath_mm = arrays['path_length_m'] * 1e3
    tau = arrays['tortuosity']
    transmitted = arrays['transmitted_in_window'].astype(bool)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    ax = axes[0, 0]
    ax.hist(fpt / 3600.0, bins=80, color='#2878b5', alpha=0.85)
    ax.axvline(observation_time_s / 3600.0, color='#c1121f', ls='--', lw=2, label='observation window')
    ax.set_xlabel('First-passage time (h)')
    ax.set_ylabel('Particles')
    ax.set_title('Eventual FPT distribution')
    ax.legend()

    ax = axes[0, 1]
    ax.hist(np.log10(lpath_mm), bins=80, color='#c65d24', alpha=0.85)
    ax.set_xlabel('log10(coarse path length / mm)')
    ax.set_ylabel('Particles')
    ax.set_title('Coarse Brownian path length')

    ax = axes[1, 0]
    ax.hist(np.log10(tau), bins=80, color='#6a994e', alpha=0.85)
    ax.set_xlabel('log10(tortuosity tau = Lpath / Lfilm)')
    ax.set_ylabel('Particles')
    ax.set_title('Tortuosity distribution')

    ax = axes[1, 1]
    labels = ['transmitted', 'not yet transmitted']
    vals = [int(transmitted.sum()), int((~transmitted).sum())]
    ax.bar(labels, vals, color=['#1b9e77', '#7570b3'])
    ax.set_ylabel('Particles')
    ax.set_title('Finite-time Ptrans count')
    fig.suptitle('Trilayer Brownian/kMC barrier metrics', fontsize=16, fontweight='bold')
    fig.savefig(out_dir / 'trilayer_metric_distributions.png', dpi=220)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--particles', type=int, default=200000)
    ap.add_argument('--seed', type=int, default=20260522)
    ap.add_argument('--out-dir', type=Path, default=Path('/home/xiao/文档/KMC/trilayer_brownian_kmc'))
    ap.add_argument('--observation-time-s', type=float, default=86400.0)
    ap.add_argument('--temp-c', type=float, default=38.0)
    ap.add_argument('--rh-high', type=float, default=0.90)
    ap.add_argument('--rh-low', type=float, default=0.0)
    ap.add_argument('--single-wvtr-mg-m2-day', type=float, default=1.26)
    ap.add_argument('--trilayer-wvtr-mg-m2-day', type=float, default=0.046)
    ap.add_argument('--d-al2o3', type=float, default=5.26e-19)
    ap.add_argument('--d-polymer', type=float, default=2.0e-13)
    ap.add_argument('--al2o3-nm', type=float, default=44.0)
    ap.add_argument('--polymer-nm', type=float, default=208.0)
    ap.add_argument('--defect-diameter-um', type=float, default=1.0)
    ap.add_argument('--defect-spacing-um', type=float, default=100.0)
    ap.add_argument('--coarse-hop-um', type=float, default=1.0)
    ap.add_argument('--gamma-shape-search', type=float, default=1.2)
    ap.add_argument('--gamma-shape-al2o3', type=float, default=2.0)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    l_al = args.al2o3_nm * 1e-9
    l_poly = args.polymer_nm * 1e-9
    l_total = 2.0 * l_al + l_poly
    defect_radius = args.defect_diameter_um * 0.5e-6
    defect_spacing = args.defect_spacing_um * 1e-6
    coarse_hop = args.coarse_hop_um * 1e-6

    # Boundary concentration is inferred from the single-layer WVTR so the
    # simulated D_eff maps back onto the reported experimental WVTR scale.
    j_single_kg_m2_s = args.single_wvtr_mg_m2_day * 1e-6 / 86400.0
    delta_c_from_single = j_single_kg_m2_s * l_al / args.d_al2o3

    d_eff_observed = (
        args.d_al2o3
        * (args.trilayer_wvtr_mg_m2_day / args.single_wvtr_mg_m2_day)
        * (l_total / l_al)
    )
    fpt_mean_target = l_total * l_total / (2.0 * d_eff_observed)
    al_mean = l_al * l_al / (2.0 * args.d_al2o3)
    search_mean = max(1e-30, fpt_mean_target - 2.0 * al_mean)

    # Theory diagnostics for defect decoupling.
    defect_area_fraction = math.pi * defect_radius * defect_radius / (defect_spacing * defect_spacing)
    mean_nearest_defect_distance = defect_spacing * (
        math.sqrt(2.0) + math.log(1.0 + math.sqrt(2.0))
    ) / 6.0
    narrow_escape_2pi = (
        defect_spacing * defect_spacing
        / (2.0 * math.pi * args.d_polymer)
        * math.log((defect_spacing / math.sqrt(math.pi)) / defect_radius)
    )
    narrow_escape_penalty = search_mean / narrow_escape_2pi

    n = args.particles
    al_scale = al_mean / args.gamma_shape_al2o3
    search_scale = search_mean / args.gamma_shape_search
    t_al_top = rng.gamma(args.gamma_shape_al2o3, al_scale, n)
    t_search = rng.gamma(args.gamma_shape_search, search_scale, n)
    t_al_bottom = rng.gamma(args.gamma_shape_al2o3, al_scale, n)
    fpt = t_al_top + t_search + t_al_bottom

    transmitted = fpt <= args.observation_time_s

    # Coarse path length at a chosen Brownian-hop resolution. In 2D Brownian
    # motion, the arc length depends on the observation scale; 1 um is chosen
    # because it matches the defect diameter scale in the schematic.
    hop_time = coarse_hop * coarse_hop / (4.0 * args.d_polymer)
    polymer_hops = np.maximum(1.0, t_search / hop_time)
    path_length_m = l_total + polymer_hops * coarse_hop
    tortuosity = path_length_m / l_total

    d_eff_sim = float(l_total * l_total / (2.0 * np.mean(fpt)))
    wvtr_sim_mg_m2_day = float(d_eff_sim * delta_c_from_single / l_total * 86400.0 * 1e6)
    wvtr_ideal_delta_c_mg_m2_day = float(
        d_eff_sim
        * (water_vapor_density_kg_m3(args.temp_c, args.rh_high) - water_vapor_density_kg_m3(args.temp_c, args.rh_low))
        / l_total
        * 86400.0
        * 1e6
    )

    metrics = {
        'model': 'coarse_grained_trilayer_brownian_kmc',
        'assumptions': {
            'description': 'Al2O3 defects are through-layer entry/exit windows; polymer interlayer decouples defects and imposes a Brownian lateral search penalty calibrated to observed trilayer WVTR.',
            'Al2O3_top_nm': args.al2o3_nm,
            'polymer_nm': args.polymer_nm,
            'Al2O3_bottom_nm': args.al2o3_nm,
            'total_thickness_nm': l_total * 1e9,
            'D_Al2O3_m2_s': args.d_al2o3,
            'D_polymer_m2_s': args.d_polymer,
            'defect_diameter_um': args.defect_diameter_um,
            'defect_spacing_um': args.defect_spacing_um,
            'coarse_path_hop_um': args.coarse_hop_um,
            'path_length_note': 'Brownian path length is scale-dependent; this model reports coarse arc length at the chosen hop scale.',
            'observation_time_s_for_Ptrans': args.observation_time_s,
            'delta_c_kg_m3_inferred_from_single_layer_WVTR': delta_c_from_single,
        },
        'calibration': {
            'single_layer_WVTR_mg_m2_day': args.single_wvtr_mg_m2_day,
            'target_trilayer_WVTR_mg_m2_day': args.trilayer_wvtr_mg_m2_day,
            'target_D_eff_m2_s': d_eff_observed,
            'target_mean_FPT_s_from_D_eff': fpt_mean_target,
            'Al2O3_defect_crossing_mean_s_each': al_mean,
            'polymer_lateral_search_mean_s_calibrated': search_mean,
            'defect_area_fraction': defect_area_fraction,
            'mean_nearest_defect_distance_um': mean_nearest_defect_distance * 1e6,
            'narrow_escape_time_2pi_s': narrow_escape_2pi,
            'calibrated_search_over_narrow_escape_2pi': narrow_escape_penalty,
        },
        'outputs': {
            'particles': n,
            'Ptrans_observation_window': float(np.mean(transmitted)),
            'transmitted_in_window': int(np.sum(transmitted)),
            'FPT_eventual_s': finite_summary(fpt),
            'FPT_transmitted_within_window_s': finite_summary(fpt[transmitted]),
            'Lpath_coarse_m': finite_summary(path_length_m),
            'Lpath_coarse_mm': finite_summary(path_length_m * 1e3),
            'tortuosity': finite_summary(tortuosity),
            'D_eff_m2_s': d_eff_sim,
            'WVTR_mg_m2_day_calibrated_delta_c': wvtr_sim_mg_m2_day,
            'WVTR_mg_m2_day_ideal_vapor_delta_c': wvtr_ideal_delta_c_mg_m2_day,
        },
    }

    with open(args.out_dir / 'trilayer_brownian_kmc_metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    rows: list[tuple[str, object]] = []
    flatten('', metrics, rows)
    with open(args.out_dir / 'trilayer_brownian_kmc_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['metric', 'value'])
        w.writerows(rows)

    arrays = {
        'first_passage_time_s': fpt,
        'transmitted_in_window': transmitted.astype(np.int8),
        'path_length_m': path_length_m,
        'tortuosity': tortuosity,
        't_al_top_s': t_al_top,
        't_polymer_search_s': t_search,
        't_al_bottom_s': t_al_bottom,
    }
    np.savez_compressed(args.out_dir / 'trilayer_brownian_kmc_arrays.npz', **arrays)
    make_plots(args.out_dir, arrays, args.observation_time_s)

    report = f"""# 三层膜 Brownian/kMC 粗粒化模型

## 文献约束后的建模思想

有机/无机多层阻隔膜的关键不是简单串联厚度，而是有机层把上下无机层缺陷解耦，
使水分子必须在有机层中横向扩散搜索下一个缺陷入口。Al2O3 层本体近似不可渗透，
渗透主要由缺陷控制。

## 本次模型

- 结构: Al2O3 / polymer / Al2O3 = {args.al2o3_nm:g} / {args.polymer_nm:g} / {args.al2o3_nm:g} nm
- D_Al2O3 = {args.d_al2o3:.3e} m2/s
- D_polymer = {args.d_polymer:.3e} m2/s
- 缺陷直径 = {args.defect_diameter_um:g} um
- 缺陷间距 = {args.defect_spacing_um:g} um
- 路径长度粗粒化步长 = {args.coarse_hop_um:g} um
- Ptrans 观测窗口 = {args.observation_time_s:g} s

## 输出

- Ptrans({args.observation_time_s:g} s) = {metrics['outputs']['Ptrans_observation_window']:.4f}
- FPT mean = {metrics['outputs']['FPT_eventual_s']['mean']:.3g} s
- FPT median = {metrics['outputs']['FPT_eventual_s']['median']:.3g} s
- Lpath mean = {metrics['outputs']['Lpath_coarse_mm']['mean']:.3g} mm
- tau mean = {metrics['outputs']['tortuosity']['mean']:.3g}
- D_eff = {metrics['outputs']['D_eff_m2_s']:.3e} m2/s
- WVTR = {metrics['outputs']['WVTR_mg_m2_day_calibrated_delta_c']:.4g} mg m^-2 day^-1

## 解释

模型被校准到图中三层膜 WVTR = {args.trilayer_wvtr_mg_m2_day:g} mg m^-2 day^-1。
由此反推出三层膜 D_eff = {d_eff_observed:.3e} m2/s，平均首次通过时间约
{fpt_mean_target:.3g} s。单个 Al2O3 缺陷层的穿越时间尺度仅约 {al_mean:.3g} s，
所以主要阻滞来自 polymer 层中的横向搜索，校准均值约 {search_mean:.3g} s。

缺陷面积分数为 {defect_area_fraction:.3e}，直接上下对齐的概率很低。
二维窄逃逸估计给出搜索时间约 {narrow_escape_2pi:.3g} s；校准搜索时间更长，
可理解为真实缺陷入口并非理想吸收边界、界面滞留/吸附、路径重复试探和局部堵塞
共同造成的额外惩罚。
"""
    (args.out_dir / 'trilayer_brownian_kmc_report.md').write_text(report, encoding='utf-8')

    print(json.dumps({
        'out_dir': str(args.out_dir),
        'Ptrans_observation_window': metrics['outputs']['Ptrans_observation_window'],
        'FPT_mean_s': metrics['outputs']['FPT_eventual_s']['mean'],
        'Lpath_mean_mm': metrics['outputs']['Lpath_coarse_mm']['mean'],
        'tau_mean': metrics['outputs']['tortuosity']['mean'],
        'D_eff_m2_s': metrics['outputs']['D_eff_m2_s'],
        'WVTR_mg_m2_day': metrics['outputs']['WVTR_mg_m2_day_calibrated_delta_c'],
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
