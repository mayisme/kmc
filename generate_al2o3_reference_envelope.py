#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from collections import deque
from pathlib import Path

import numpy as np

DX_NM = 0.5
THICKNESS_NM = 10.0
WIDTH_NM = 320.0
TEMP_C = 38.0
RH_HIGH = 0.90
LITERATURE_SORBED_C1_KG_M3 = 6.77


def saturation_vapor_pressure_pa(temp_c: float) -> float:
    return 611.21 * math.exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))


def water_vapor_density_kg_m3(temp_c: float, rh: float) -> float:
    p = saturation_vapor_pressure_pa(temp_c) * rh
    return p * 0.01801528 / (8.314462618 * (temp_c + 273.15))


def deff_from_wvtr(wvtr_g_m2_day: float, film_thickness_nm: float = THICKNESS_NM,
                   delta_c_kg_m3: float = LITERATURE_SORBED_C1_KG_M3) -> float:
    return wvtr_g_m2_day * (film_thickness_nm * 1e-9) / (86400.0 * 1000.0 * delta_c_kg_m3)


def wvtr_ideal_from_deff(deff_m2_s: float, film_thickness_nm: float = THICKNESS_NM) -> float:
    delta_c = water_vapor_density_kg_m3(TEMP_C, RH_HIGH)
    return deff_m2_s * delta_c / (film_thickness_nm * 1e-9) * 86400.0 * 1000.0


def draw_disk(grid: np.ndarray, y0: float, x0: float, radius_cells: float, value: int) -> None:
    h, w = grid.shape
    r2 = radius_cells * radius_cells
    ymin = max(0, int(math.floor(y0 - radius_cells - 1)))
    ymax = min(h, int(math.ceil(y0 + radius_cells + 2)))
    xmin = int(math.floor(x0 - radius_cells - 1))
    xmax = int(math.ceil(x0 + radius_cells + 2))
    for y in range(ymin, ymax):
        for xx in range(xmin, xmax):
            x = xx % w
            dx = min(abs(xx - x0), w - abs(xx - x0))
            dy = y - y0
            if dx * dx + dy * dy <= r2:
                grid[y, x] = value


def build_staggered_nonthrough_defects(seed: int = 20260521):
    rng = np.random.default_rng(seed)
    width = int(round(WIDTH_NM / DX_NM))
    height = int(round(THICKNESS_NM / DX_NM))
    grid = np.zeros((height, width), dtype=np.int32)
    top_mask = np.zeros(width, dtype=np.int32)
    defects = []

    layer_slices = [(0, 4), (5, 9), (11, 14), (16, 19)]
    base_centers_nm = []
    x_nm = rng.uniform(0, 70.0)
    while x_nm < WIDTH_NM:
        base_centers_nm.append(float(x_nm))
        x_nm += rng.uniform(55.0, 85.0)

    for li, (y0, y1) in enumerate(layer_slices):
        offset_nm = [0.0, 28.0, -24.0, 36.0][li]
        for c_nm in base_centers_nm:
            if rng.random() < 0.18:
                continue
            cx = ((c_nm + offset_nm + rng.normal(0.0, 2.0)) % WIDTH_NM) / DX_NM
            d_nm = rng.uniform(1.0, 2.0)
            radius_cells = max(1.0, (d_nm / 2.0) / DX_NM)
            for yy in range(y0, y1 + 1):
                draw_disk(grid, yy, cx + rng.normal(0.0, 0.20), radius_cells, 1)
            if y0 == 0:
                cx_int = int(round(cx)) % width
                r_int = max(1, int(round(radius_cells)))
                for dx_i in range(-r_int, r_int + 1):
                    top_mask[(cx_int + dx_i) % width] = 1
            defects.append({
                'layer_index': li,
                'y_start_cell': y0,
                'y_end_cell': y1,
                'center_nm': float((cx * DX_NM) % WIDTH_NM),
                'diameter_nm': float(d_nm),
            })

    return grid, top_mask, defects


def connectivity(mask: np.ndarray):
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    for y in range(h):
        for x in range(w):
            if not mask[y, x] or seen[y, x]:
                continue
            q = deque([(y, x)])
            seen[y, x] = True
            count = 0
            top = False
            bottom = False
            min_y = h
            max_y = -1
            xs = set()
            while q:
                yy, xx = q.popleft()
                count += 1
                top |= yy == 0
                bottom |= yy == h - 1
                min_y = min(min_y, yy)
                max_y = max(max_y, yy)
                xs.add(xx)
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny = yy + dy
                    nx = (xx + dx) % w
                    if 0 <= ny < h and mask[ny, nx] and not seen[ny, nx]:
                        seen[ny, nx] = True
                        q.append((ny, nx))
            comps.append({
                'cells': count,
                'touches_top': top,
                'touches_bottom': bottom,
                'min_y': min_y,
                'max_y': max_y,
                'x_span_cells': len(xs),
            })
    percolating = [c for c in comps if c['touches_top'] and c['touches_bottom']]
    return comps, percolating


def write_json_csv_npz(out_dir: Path, metrics: dict, arrays: dict[str, np.ndarray]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / 'kmc_transport_metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    def flatten(prefix: str, obj, rows: list[tuple[str, object]]) -> None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                flatten(f'{prefix}.{k}' if prefix else k, v, rows)
        elif isinstance(obj, list):
            rows.append((prefix, json.dumps(obj, ensure_ascii=False)))
        else:
            rows.append((prefix, obj))

    rows: list[tuple[str, object]] = []
    flatten('', metrics, rows)
    with open(out_dir / 'kmc_transport_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['metric', 'value'])
        w.writerows(rows)
    if arrays:
        np.savez_compressed(out_dir / 'kmc_transport_arrays.npz', **arrays)


def finite_empty_summary():
    return {
        'count': 0, 'mean': None, 'median': None, 'std': None,
        'p01': None, 'p05': None, 'p10': None, 'p25': None,
        'p75': None, 'p90': None, 'p95': None, 'p99': None,
        'min': None, 'max': None,
    }


def metrics_template(source: str, model: str, ptrans: float, deff: float | None,
                     wvtr_lit: float | None, wvtr_ideal: float | None,
                     assumptions_extra: dict):
    assumptions = {
        'dx_nm': DX_NM,
        'film_thickness_nm': THICKNESS_NM,
        'temperature_c': TEMP_C,
        'rh_high': RH_HIGH,
        'rh_low': 0.0,
        'wvtr_formula': 'WVTR = D_eff * delta_C / L * 86400 * 1000, units g m^-2 day^-1',
        'literature_sorbed_delta_c_kg_m3_from_C1_0.00677_g_cm3': LITERATURE_SORBED_C1_KG_M3,
    }
    assumptions.update(assumptions_extra)
    return {
        'source_raw_npz': source,
        'model': model,
        'assumptions': assumptions,
        'Ptrans': ptrans,
        'particles': None,
        'transmitted': None,
        'non_transmitted': None,
        'FPT_success_s': finite_empty_summary(),
        'path_length_success_mm': finite_empty_summary(),
        'tortuosity_success': finite_empty_summary(),
        'D_eff_m2_s': {
            'success_mean_FPT': deff,
            'success_median_FPT': None,
        },
        'WVTR_g_m2_day': {
            'using_ideal_vapor_delta_c': wvtr_ideal,
            'using_user_delta_c': None,
            'using_literature_sorbed_C1_6p77_kg_m3': wvtr_lit,
            'relative_WVTR_factor_equals_relative_D_eff': True,
        },
    }


def make_defect_plot(out_dir: Path, grid: np.ndarray) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(14, 4), constrained_layout=True)
    ax.imshow(grid == 1, aspect='auto', interpolation='nearest', cmap='gray_r')
    ax.set_title('10 nm Al2O3 staggered non-through defect reference')
    ax.set_xlabel('x cell')
    ax.set_ylabel('z cell')
    fig.savefig(out_dir / 'staggered_nonthrough_defect_mask.png', dpi=220)
    plt.close(fig)


def main() -> None:
    root = Path('/home/xiao/文档/KMC')
    out_root = root / 'al2o3_reference_envelope'
    out_root.mkdir(parents=True, exist_ok=True)

    grid, top_mask, defects = build_staggered_nonthrough_defects()
    comps, percolating = connectivity(grid == 1)
    nonthrough_dir = out_root / 'al2o3_10nm_staggered_nonthrough'
    nonthrough_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(nonthrough_dir / 'al2o3_10nm_staggered_nonthrough_raw.npz',
                        grid=grid, pinhole_mask=top_mask)
    metrics = metrics_template(
        source=str(nonthrough_dir / 'al2o3_10nm_staggered_nonthrough_raw.npz'),
        model='single_al2o3_10nm_staggered_nonthrough_defects_connectivity_limited',
        ptrans=0.0,
        deff=0.0,
        wvtr_lit=0.0,
        wvtr_ideal=0.0,
        assumptions_extra={
            'transport_assumption': 'Dense Al2O3 matrix is impermeable and staggered defects do not form a top-to-bottom connected path; transport is connectivity-limited to zero in this idealized model.',
            'defect_cell_fraction': float((grid == 1).mean()),
            'top_mask_fraction': float(top_mask.mean()),
            'connected_defect_components': len(comps),
            'top_to_bottom_percolating_defect_components': len(percolating),
            'defects': defects,
        },
    )
    write_json_csv_npz(nonthrough_dir / 'transport_metrics', metrics, {})
    make_defect_plot(nonthrough_dir, grid)

    literature_cases = [
        {
            'id': 'dense_uv_ald_10nm_pet_mocon',
            'label': '10 nm UV-ALD Al2O3/PET, MOCON literature-calibrated',
            'wvtr_lit': 1.42e-3,
            'source': 'RSC Advances 2017, UV-enhanced ALD Al2O3 moisture barrier; reported 10 nm UV-ALD Al2O3 WVTR around 1.42e-3 g m^-2 day^-1.',
        },
        {
            'id': 'thermal_ald_10nm_pet_mocon',
            'label': '10 nm thermal ALD Al2O3/PET, poorer-process literature-calibrated',
            'wvtr_lit': 6.11e-1,
            'source': 'RSC Advances 2017, thermal ALD 10 nm Al2O3 comparison; reported WVTR around 6.11e-1 g m^-2 day^-1.',
        },
        {
            'id': 'pa_ald_20nm_pen_context',
            'label': '20 nm PA-ALD Al2O3/PEN context value',
            'wvtr_lit': 5.0e-3,
            'source': 'Plasma-assisted ALD Al2O3/PEN moisture-barrier study; reported WVTR around 5e-3 g m^-2 day^-1 for 20 nm Al2O3.',
        },
    ]

    for case in literature_cases:
        deff = deff_from_wvtr(case['wvtr_lit'])
        metrics = metrics_template(
            source='literature_calibrated_no_raw_kmc',
            model=case['id'],
            ptrans=float('nan'),
            deff=deff,
            wvtr_lit=case['wvtr_lit'],
            wvtr_ideal=wvtr_ideal_from_deff(deff),
            assumptions_extra={
                'transport_assumption': 'Literature-calibrated continuum reference. Ptrans/FPT/Lpath/tau are not defined without a microscopic trajectory model.',
                'label': case['label'],
                'literature_note': case['source'],
                'calibration_delta_c_kg_m3': LITERATURE_SORBED_C1_KG_M3,
            },
        )
        write_json_csv_npz(out_root / case['id'] / 'transport_metrics', metrics, {})

    summary = {
        'generated': str(out_root),
        'references': [case['id'] for case in literature_cases],
        'nonthrough_reference': 'al2o3_10nm_staggered_nonthrough',
        'note': 'Use these references as an uncertainty envelope, not as one-to-one KMC trajectory equivalents.',
    }
    (out_root / 'reference_envelope_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
