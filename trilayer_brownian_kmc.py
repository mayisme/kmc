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


def make_sorption_plot(out_dir: Path, fpt_steady: np.ndarray, fpt_effective: np.ndarray,
                       sorption: dict, observation_time_s: float) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    tau = sorption['t_saturation_s']
    bt_days = sorption['breakthrough_time_days']
    t50d = sorption['t50_saturation_days']
    t90d = sorption['t90_saturation_days']
    t99d = sorption['t99_saturation_days']
    obs_days = observation_time_s / 86400.0

    # Time axis out to ~1.2 * (t99 or longest FPT), in days.
    horizon_days = max(t99d * 1.15, float(np.percentile(fpt_effective, 99)) / 86400.0)
    if not np.isfinite(horizon_days) or horizon_days <= 0:
        horizon_days = 1.0
    t_days = np.linspace(0.0, horizon_days, 600)
    t_s = t_days * 86400.0

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)

    # Panel 1: one-sided sorption fill toward FULL saturation, sat(t) = 1 - exp(-t/t_sat).
    ax = axes[0]
    if np.isfinite(tau) and tau > 0:
        sat = 1.0 - np.exp(-t_s / tau)
    else:
        sat = np.zeros_like(t_s)
    ax.plot(t_days, sat * 100.0, color='#2878b5', lw=2.4, label='polymer water saturation')
    for td, lbl, col in [(t50d, 't50', '#6a994e'), (t90d, 't90', '#e09f3e'), (t99d, 't99', '#c1121f')]:
        if np.isfinite(td) and td <= horizon_days:
            ax.axvline(td, color=col, ls=':', lw=1.6)
            ax.text(td, 5, f'{lbl}\n{td:.1f} d', color=col, fontsize=9, ha='center')
    ax.axvline(bt_days, color='#9d4edd', ls='--', lw=2.0,
               label=f'saturation t_sat = M_sat/J_in = {bt_days:.1f} d')
    if obs_days <= horizon_days:
        ax.axvline(obs_days, color='#333333', ls='-', lw=1.4, alpha=0.7,
                   label=f'observation window = {obs_days:.2f} d')
    ax.set_xlabel('time (days)')
    ax.set_ylabel('polymer saturation (%)')
    ax.set_title('Water-vapour sorption fill (throttled by top Al2O3)')
    ax.set_ylim(0, 105)
    ax.legend(loc='lower right', fontsize=9)

    # Panel 2: cumulative transmission Ptrans(t), steady vs sorption-gated.
    ax = axes[1]
    fs = np.sort(fpt_steady[np.isfinite(fpt_steady)])
    fe = np.sort(fpt_effective[np.isfinite(fpt_effective)])
    p_steady = np.searchsorted(fs, t_s, side='right') / max(1, fs.size)
    p_eff = np.searchsorted(fe, t_s, side='right') / max(1, fe.size)
    ax.plot(t_days, p_steady * 100.0, color='#c65d24', lw=2.2, label='no sorption (steady only)')
    ax.plot(t_days, p_eff * 100.0, color='#1b9e77', lw=2.4, label='with sorption breakthrough')
    ax.axvline(bt_days, color='#9d4edd', ls='--', lw=2.0, label=f'breakthrough = {bt_days:.1f} d')
    if obs_days <= horizon_days:
        ax.axvline(obs_days, color='#333333', ls='-', lw=1.4, alpha=0.7,
                   label=f'observation window = {obs_days:.2f} d')
    ax.set_xlabel('time (days)')
    ax.set_ylabel('cumulative transmitted fraction Ptrans(t) (%)')
    ax.set_title('Breakthrough delay of finite-time transmission')
    ax.legend(loc='lower right', fontsize=9)

    fig.suptitle('Trilayer water-vapour sorption breakthrough (saturate-then-cross)',
                 fontsize=15, fontweight='bold')
    fig.savefig(out_dir / 'trilayer_sorption_breakthrough.png', dpi=220)
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
    ap.add_argument('--solubility-mol-m3-pa', type=float, default=1.0,
                    help='Polymer water solubility coefficient S (Henry law): C_eq = S * p_upstream [mol m^-3 Pa^-1]')
    ap.add_argument('--sorption-lag', action=argparse.BooleanOptionalAction, default=True,
                    help='Include the solution-diffusion sorption breakthrough (RC / Daynes-Barrer time-lag) gate (Option A)')
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

    # --- Solution-diffusion sorption reservoir (Option A: RC / Daynes-Barrer time-lag) ---
    # The polymer interlayer is a Henry-law sorption "capacitor" that is charged through
    # the upstream Al2O3 (a high resistance) and discharged through the downstream Al2O3.
    # Until this reservoir is charged, downstream breakthrough is gated, producing a
    # transient induction time. This does NOT change the steady-state, defect-limited WVTR;
    # it only delays the transient (Ptrans(t)) — exactly why high-solubility interlayers
    # extend the breakthrough/lag time without lowering steady-state permeation.
    p_sat = saturation_vapor_pressure_pa(args.temp_c)
    p_up = p_sat * args.rh_high
    p_down = p_sat * args.rh_low
    delta_p = p_up - p_down
    c_eq = args.solubility_mol_m3_pa * p_up                 # mol/m^3, equilibrium dissolved conc
    areal_capacitance = args.solubility_mol_m3_pa * l_poly  # mol/(m^2.Pa), areal sorption capacitance
    m_sat_mol_m2 = c_eq * l_poly                            # mol/m^2 to fully saturate the layer
    m_sat_mg_m2 = m_sat_mol_m2 * 0.01801528 * 1e6           # mg/m^2
    # Single-Al2O3 areal permeance W inferred from the measured single-layer WVTR.
    j_single_mol_m2_s = j_single_kg_m2_s / 0.01801528
    w_al2o3 = j_single_mol_m2_s / p_up if p_up > 0.0 else 0.0   # mol/(m^2.s.Pa), top-Al2O3 permeance
    # Water-vapour permeation, sequential gate: a molecule enters through a top-Al2O3 defect,
    # then the polymer must absorb water vapour up to its Henry saturation (C_eq) BEFORE
    # crossing the bottom Al2O3 begins. The fill is supplied (throttled) by the sparse
    # top-Al2O3 water-vapour influx; the bottom is gated shut during the fill, so the layer
    # charges ONE-SIDEDLY toward FULL saturation (no bottom "discharge").
    # Initial influx J_in = W*p1 = single-layer molar flux; saturation time t_sat = M_sat/J_in = C/W.
    t_sat_s = areal_capacitance / w_al2o3 if w_al2o3 > 0.0 else float('inf')   # = M_sat / J_in
    # One-sided sorption fill toward full C_eq: sat(t) = 1 - exp(-t/t_sat).
    t50_s = t_sat_s * math.log(2.0)
    t90_s = t_sat_s * math.log(10.0)
    t99_s = t_sat_s * math.log(100.0)
    # Pure-diffusion (Daynes-Barrer) lag of the polymer alone; negligible here, kept for context.
    diffusive_time_lag_s = l_poly * l_poly / (6.0 * args.d_polymer)
    breakthrough_time_s = t_sat_s   # time to saturate the polymer = water-vapour breakthrough lag
    sorption_lag_s = breakthrough_time_s if args.sorption_lag else 0.0

    n = args.particles
    al_scale = al_mean / args.gamma_shape_al2o3
    search_scale = search_mean / args.gamma_shape_search
    t_al_top = rng.gamma(args.gamma_shape_al2o3, al_scale, n)
    t_search = rng.gamma(args.gamma_shape_search, search_scale, n)
    t_al_bottom = rng.gamma(args.gamma_shape_al2o3, al_scale, n)
    fpt_steady = t_al_top + t_search + t_al_bottom

    # Option A: additive series induction — the sorption reservoir first charges (breakthrough
    # lag), then steady defect-limited crossing proceeds. Steady-state metrics use fpt_steady;
    # the transient (Ptrans) uses the sorption-gated effective fpt.
    fpt = fpt_steady + sorption_lag_s
    transmitted = fpt <= args.observation_time_s
    transmitted_no_sorption = fpt_steady <= args.observation_time_s

    # Coarse path length at a chosen Brownian-hop resolution. In 2D Brownian
    # motion, the arc length depends on the observation scale; 1 um is chosen
    # because it matches the defect diameter scale in the schematic.
    hop_time = coarse_hop * coarse_hop / (4.0 * args.d_polymer)
    polymer_hops = np.maximum(1.0, t_search / hop_time)
    path_length_m = l_total + polymer_hops * coarse_hop
    tortuosity = path_length_m / l_total

    d_eff_sim = float(l_total * l_total / (2.0 * np.mean(fpt_steady)))
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
        'sorption_breakthrough': {
            'model': 'water_vapor_permeation_one_sided_sorption_fill_to_saturation',
            'enabled': bool(args.sorption_lag),
            'solubility_S_mol_m3_pa': args.solubility_mol_m3_pa,
            'temp_c': args.temp_c,
            'rh_high': args.rh_high,
            'rh_low': args.rh_low,
            'p_sat_pa': p_sat,
            'p_upstream_pa': p_up,
            'p_downstream_pa': p_down,
            'delta_p_pa': delta_p,
            'C_eq_mol_m3': c_eq,
            'areal_capacitance_mol_m2_pa': areal_capacitance,
            'M_saturation_mol_m2': m_sat_mol_m2,
            'M_saturation_mg_m2': m_sat_mg_m2,
            'top_Al2O3_water_vapor_influx_mol_m2_s': j_single_mol_m2_s,
            'top_Al2O3_permeance_mol_m2_s_pa': w_al2o3,
            't_saturation_s': t_sat_s,
            't_saturation_days': t_sat_s / 86400.0,
            'breakthrough_time_s': breakthrough_time_s,
            'breakthrough_time_days': breakthrough_time_s / 86400.0,
            't50_saturation_s': t50_s,
            't50_saturation_days': t50_s / 86400.0,
            't90_saturation_s': t90_s,
            't90_saturation_days': t90_s / 86400.0,
            't99_saturation_s': t99_s,
            't99_saturation_days': t99_s / 86400.0,
            'diffusive_time_lag_polymer_s': diffusive_time_lag_s,
            'sorption_lag_applied_s': sorption_lag_s,
            'note': 'Sequential water-vapour gate: the polymer must absorb water to full Henry saturation (C_eq = S*p1), throttled by the sparse top-Al2O3 influx, BEFORE bottom-Al2O3 crossing begins. t_saturation = M_sat/J_in = C/W (one-sided fill, no bottom discharge). This delays the transient (Ptrans) only; the steady-state defect-limited WVTR/D_eff are computed from the steady FPT and are unchanged.',
        },
        'outputs': {
            'particles': n,
            'Ptrans_observation_window': float(np.mean(transmitted)),
            'transmitted_in_window': int(np.sum(transmitted)),
            'Ptrans_observation_window_no_sorption': float(np.mean(transmitted_no_sorption)),
            'transmitted_in_window_no_sorption': int(np.sum(transmitted_no_sorption)),
            'FPT_eventual_s': finite_summary(fpt),
            'FPT_steady_no_sorption_s': finite_summary(fpt_steady),
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
        'first_passage_time_steady_s': fpt_steady,
        'transmitted_in_window': transmitted.astype(np.int8),
        'transmitted_in_window_no_sorption': transmitted_no_sorption.astype(np.int8),
        'path_length_m': path_length_m,
        'tortuosity': tortuosity,
        't_al_top_s': t_al_top,
        't_polymer_search_s': t_search,
        't_al_bottom_s': t_al_bottom,
    }
    np.savez_compressed(args.out_dir / 'trilayer_brownian_kmc_arrays.npz', **arrays)
    make_plots(args.out_dir, arrays, args.observation_time_s)
    sorption_info = metrics['sorption_breakthrough']
    make_sorption_plot(args.out_dir, fpt_steady, fpt, sorption_info, args.observation_time_s)

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

## 输出（稳态，缺陷控制；不受吸附闸门影响）

- Ptrans({args.observation_time_s:g} s) = {metrics['outputs']['Ptrans_observation_window']:.4f}  （含吸附击穿闸门）
- Ptrans({args.observation_time_s:g} s) = {metrics['outputs']['Ptrans_observation_window_no_sorption']:.4f}  （旧模型，无吸附）
- FPT(steady) mean = {metrics['outputs']['FPT_steady_no_sorption_s']['mean']:.3g} s
- FPT(effective, 含击穿) mean = {metrics['outputs']['FPT_eventual_s']['mean']:.3g} s
- Lpath mean = {metrics['outputs']['Lpath_coarse_mm']['mean']:.3g} mm
- tau mean = {metrics['outputs']['tortuosity']['mean']:.3g}
- D_eff(steady) = {metrics['outputs']['D_eff_m2_s']:.3e} m2/s
- WVTR(steady) = {metrics['outputs']['WVTR_mg_m2_day_calibrated_delta_c']:.4g} mg m^-2 day^-1

## 溶解度驱动的吸水饱和击穿（水蒸气渗透：先饱和后穿透）

水蒸气分子穿过上层 Al2O3 缺陷后进入高分子层。按本项目设定的顺序机理：**高分子层必须先把
水吸到 Henry 饱和（C_eq = S·p1），其供水受上层 Al2O3 稀疏缺陷的水蒸气进入通量节流；在达到
饱和之前，下层 Al2O3 的穿透尚未开始**（单向充填至满饱和，无"下层放电"）。**该机制只延迟瞬态
（Ptrans(t)），不改变由缺陷控制的稳态 WVTR/D_eff**——这正是高溶解度有机夹层能拉长击穿/滞后
时间却不降低稳态渗透的原因。

- 溶解度系数 S = {sorption_info['solubility_S_mol_m3_pa']:g} mol m^-3 Pa^-1
- 工况: T = {sorption_info['temp_c']:g} °C, 上游 RH = {sorption_info['rh_high']:g}, 下游 RH = {sorption_info['rh_low']:g}
- 饱和蒸气压 p_sat = {sorption_info['p_sat_pa']:.4g} Pa；上游水分压 p1 = {sorption_info['p_upstream_pa']:.4g} Pa
- 高分子内平衡溶解浓度 C_eq = S·p1 = {sorption_info['C_eq_mol_m3']:.4g} mol m^-3
- 单位面积饱和储水量 M_sat = C_eq·L_poly = {sorption_info['M_saturation_mg_m2']:.4g} mg m^-2
- 上层 Al2O3 水蒸气进入通量 J_in = {sorption_info['top_Al2O3_water_vapor_influx_mol_m2_s']:.4g} mol m^-2 s^-1
- **吸水饱和时间 t_sat = M_sat/J_in = C/W = {sorption_info['t_saturation_s']:.4g} s = {sorption_info['t_saturation_days']:.3g} 天**
- 饱和充填特征时间: t50 = {sorption_info['t50_saturation_days']:.3g} 天, t90 = {sorption_info['t90_saturation_days']:.3g} 天, t99 = {sorption_info['t99_saturation_days']:.3g} 天
- 纯扩散时间滞后 L_poly^2/6D = {sorption_info['diffusive_time_lag_polymer_s']:.3g} s（可忽略，膜太薄）

**物理结论**: 真正控制饱和的是上层缺陷稀疏的水蒸气进入节流，而非高分子内扩散。
因此吸水饱和击穿约 {sorption_info['breakthrough_time_days']:.3g} 天，远大于观测窗 {args.observation_time_s/86400:g} 天，
使 Ptrans(观测窗) 从 {metrics['outputs']['Ptrans_observation_window_no_sorption']:.4f}（无吸附）降至
{metrics['outputs']['Ptrans_observation_window']:.4f}（含吸附）；而稳态 WVTR 几乎不变。
注：单向充填随膜内浓度上升、驱动力 (p1-p_poly) 减小而趋缓，sat(t)=1-exp(-t/t_sat) 渐近至满饱和；
"达到饱和才穿过"是其顺序门近似，上面给出 t50/t90/t99 供按阈值取用。

## 解释（稳态阻滞机理）

模型被校准到图中三层膜 WVTR = {args.trilayer_wvtr_mg_m2_day:g} mg m^-2 day^-1。
由此反推出三层膜 D_eff = {d_eff_observed:.3e} m2/s，平均首次通过时间约
{fpt_mean_target:.3g} s。单个 Al2O3 缺陷层的穿越时间尺度仅约 {al_mean:.3g} s，
所以主要阻滞来自 polymer 层中的横向搜索，校准均值约 {search_mean:.3g} s。

缺陷面积分数为 {defect_area_fraction:.3e}，直接上下对齐的概率很低。
二维窄逃逸估计给出搜索时间约 {narrow_escape_2pi:.3g} s；校准搜索时间更长，
可理解为真实缺陷入口并非理想吸收边界、界面滞留/吸附、路径重复试探和局部堵塞
共同造成的额外惩罚。

## 可视化

### 吸水饱和击穿（先饱和后穿透）

![吸水饱和与击穿延迟](trilayer_sorption_breakthrough.png)

左：溶解度驱动的高分子单向吸水充填曲线 sat(t)=1-exp(-t/t_sat) 及 t50/t90/t99；
右：有限时间累计透过率 Ptrans(t)——含吸水饱和击穿 vs 旧（仅稳态）的对比。

### 指标分布

![指标分布](trilayer_metric_distributions.png)

FPT、粗粒化路径长度、迂曲度分布与有限时间透过计数。
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
