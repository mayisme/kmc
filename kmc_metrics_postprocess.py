#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def saturation_vapor_pressure_pa(temp_c: float) -> float:
    """Buck equation over water, valid near room temperature."""
    return 611.21 * math.exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))


def water_vapor_density_kg_m3(temp_c: float, rh: float) -> float:
    # Ideal gas density rho = p M / R T
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


def make_plots(out_dir: Path, metrics: dict, arrays: dict[str, np.ndarray]) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ok = arrays['transmitted'].astype(bool)
    fpt = arrays['first_times'][ok]
    path_nm = arrays['path_length_nm_success']
    tau = arrays['tortuosity_success']
    final_depth_nm = arrays['final_depth_nm']
    non_final = final_depth_nm[~ok]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    ax = axes[0, 0]
    if fpt.size:
        ax.hist(fpt / 60.0, bins=70, color='#2878b5', alpha=0.85)
        ax.axvline(np.median(fpt) / 60.0, color='black', lw=2, ls='--', label='median')
        ax.set_xlabel('First passage time (min)')
        ax.set_ylabel('Successful particles')
        ax.set_title('FPT distribution')
        ax.legend()

    ax = axes[0, 1]
    if path_nm.size:
        ax.hist(np.log10(path_nm), bins=70, color='#c65d24', alpha=0.85)
        ax.axvline(np.log10(np.median(path_nm)), color='black', lw=2, ls='--', label='median')
        ax.set_xlabel('log10(path length / nm)')
        ax.set_ylabel('Successful particles')
        ax.set_title('Water molecule path length distribution')
        ax.legend()

    ax = axes[1, 0]
    if tau.size:
        ax.hist(np.log10(tau), bins=70, color='#6a994e', alpha=0.85)
        ax.axvline(np.log10(np.median(tau)), color='black', lw=2, ls='--', label='median')
        ax.set_xlabel('log10(tortuosity tau = Lpath / Lfilm)')
        ax.set_ylabel('Successful particles')
        ax.set_title('Tortuosity distribution')
        ax.legend()

    ax = axes[1, 1]
    if non_final.size:
        ax.hist(non_final, bins=np.arange(-0.25, np.max(final_depth_nm) + 0.75, 0.5), color='#7b2cbf', alpha=0.85)
        ax.set_xlabel('Final depth of non-transmitted particles (nm)')
        ax.set_ylabel('Particles')
        ax.set_title('Trapped / censored final depth distribution')

    fig.suptitle('kMC water vapor barrier metrics: Ptrans, FPT, Lpath, tau', fontsize=16, fontweight='bold')
    fig.savefig(out_dir / 'kmc_metric_distributions.png', dpi=220)
    plt.close(fig)

    if fpt.size and path_nm.size:
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        sample_n = min(12000, fpt.size)
        idx = np.linspace(0, fpt.size - 1, sample_n).astype(int)
        ax.scatter(path_nm[idx] / 1e6, fpt[idx] / 60.0, s=4, alpha=0.25, color='#264653')
        ax.set_xlabel('Path length (mm)')
        ax.set_ylabel('FPT (min)')
        ax.set_title('FPT vs path length, successful particles')
        fig.savefig(out_dir / 'kmc_fpt_vs_path_length.png', dpi=220)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('raw_npz', type=Path)
    ap.add_argument('--out-dir', type=Path, default=None)
    ap.add_argument('--dx-nm', type=float, default=0.5)
    ap.add_argument('--film-thickness-nm', type=float, default=10.5)
    ap.add_argument('--temp-c', type=float, default=38.0)
    ap.add_argument('--rh-high', type=float, default=0.90)
    ap.add_argument('--rh-low', type=float, default=0.0)
    ap.add_argument('--sorbed-delta-c-kg-m3', type=float, default=None,
                    help='Optional calibrated sorbed water concentration difference inside film. If omitted, ideal vapor gas density is used.')
    args = ap.parse_args()

    raw = args.raw_npz.resolve()
    out_dir = args.out_dir or raw.parent / 'metrics_postprocess'
    out_dir.mkdir(parents=True, exist_ok=True)

    with np.load(raw) as data:
        first_times = data['first_times'].astype(np.float64)
        steps = data['steps'].astype(np.int64)
        transmitted = data['transmitted'].astype(np.int8)
        final_y = data['final_y'].astype(np.int64)

    ok = transmitted.astype(bool)
    n = int(transmitted.size)
    n_ok = int(ok.sum())
    path_length_nm = steps.astype(np.float64) * args.dx_nm
    tortuosity = path_length_nm / args.film_thickness_nm
    final_depth_nm = final_y.astype(np.float64) * args.dx_nm

    fpt_success = first_times[ok]
    path_success = path_length_nm[ok]
    tau_success = tortuosity[ok]

    # D_eff estimates. Main value uses successful first-passage mean to match existing model output.
    L_m = args.film_thickness_nm * 1e-9
    D_eff_success_mean = None
    D_eff_success_median = None
    if n_ok:
        D_eff_success_mean = float(L_m * L_m / (2.0 * np.mean(fpt_success)))
        D_eff_success_median = float(L_m * L_m / (2.0 * np.median(fpt_success)))

    delta_c_vapor = water_vapor_density_kg_m3(args.temp_c, args.rh_high) - water_vapor_density_kg_m3(args.temp_c, args.rh_low)
    delta_c_used = args.sorbed_delta_c_kg_m3 if args.sorbed_delta_c_kg_m3 is not None else delta_c_vapor
    wvtr = None
    if D_eff_success_mean is not None:
        wvtr = float(D_eff_success_mean * delta_c_used / L_m * 86400.0 * 1000.0)

    # Literature paper used C1=0.00677 g/cm3 = 6.77 kg/m3 as fitted/sorbed concentration for PET at 38C/90RH.
    literature_delta_c = 6.77
    wvtr_literature_c = None
    if D_eff_success_mean is not None:
        wvtr_literature_c = float(D_eff_success_mean * literature_delta_c / L_m * 86400.0 * 1000.0)

    metrics = {
        'source_raw_npz': str(raw),
        'assumptions': {
            'dx_nm': args.dx_nm,
            'film_thickness_nm': args.film_thickness_nm,
            'temperature_c': args.temp_c,
            'rh_high': args.rh_high,
            'rh_low': args.rh_low,
            'wvtr_formula': 'WVTR = D_eff * delta_C / L * 86400 * 1000, units g m^-2 day^-1',
            'delta_c_vapor_kg_m3_ideal_gas': delta_c_vapor,
            'delta_c_used_kg_m3': delta_c_used,
            'literature_sorbed_delta_c_kg_m3_from_C1_0.00677_g_cm3': literature_delta_c,
            'path_length_definition': 'Lpath = path_steps * dx_nm; includes surface and bulk lattice hops.',
            'tortuosity_definition': 'tau = Lpath / film_thickness_nm.',
            'censoring_note': 'FPT, Lpath, tau summaries are for transmitted particles unless explicitly named all_particles.',
        },
        'Ptrans': float(n_ok / n) if n else 0.0,
        'particles': n,
        'transmitted': n_ok,
        'non_transmitted': int(n - n_ok),
        'FPT_success_s': finite_summary(fpt_success),
        'FPT_all_recorded_s': finite_summary(first_times),
        'path_steps_success': finite_summary(steps[ok]),
        'path_length_success_nm': finite_summary(path_success),
        'path_length_success_um': finite_summary(path_success / 1e3),
        'path_length_success_mm': finite_summary(path_success / 1e6),
        'tortuosity_success': finite_summary(tau_success),
        'final_depth_non_transmitted_nm': finite_summary(final_depth_nm[~ok]),
        'D_eff_m2_s': {
            'success_mean_FPT': D_eff_success_mean,
            'success_median_FPT': D_eff_success_median,
        },
        'WVTR_g_m2_day': {
            'using_ideal_vapor_delta_c': wvtr if args.sorbed_delta_c_kg_m3 is None else None,
            'using_user_delta_c': wvtr if args.sorbed_delta_c_kg_m3 is not None else None,
            'using_literature_sorbed_C1_6p77_kg_m3': wvtr_literature_c,
            'relative_WVTR_factor_equals_relative_D_eff': True,
        },
    }

    with open(out_dir / 'kmc_transport_metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    rows: list[tuple[str, object]] = []
    flatten('', metrics, rows)
    with open(out_dir / 'kmc_transport_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['metric', 'value'])
        w.writerows(rows)

    np.savez_compressed(
        out_dir / 'kmc_transport_arrays.npz',
        transmitted=transmitted,
        first_times_s=first_times,
        path_steps=steps,
        path_length_nm=path_length_nm,
        tortuosity=tortuosity,
        final_depth_nm=final_depth_nm,
        success_mask=ok,
        path_length_nm_success=path_success,
        tortuosity_success=tau_success,
        fpt_success_s=fpt_success,
    )

    make_plots(out_dir, metrics, {
        'transmitted': transmitted,
        'first_times': first_times,
        'path_length_nm_success': path_success,
        'tortuosity_success': tau_success,
        'final_depth_nm': final_depth_nm,
    })

    print(json.dumps({
        'out_dir': str(out_dir),
        'Ptrans': metrics['Ptrans'],
        'FPT_success_s': metrics['FPT_success_s'],
        'path_length_success_mm': metrics['path_length_success_mm'],
        'tortuosity_success': metrics['tortuosity_success'],
        'D_eff_m2_s': metrics['D_eff_m2_s'],
        'WVTR_g_m2_day': metrics['WVTR_g_m2_day'],
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
