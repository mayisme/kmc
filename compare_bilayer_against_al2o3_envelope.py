#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path


def load(path: Path) -> dict:
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def metric(m: dict, dotted: str):
    cur = m
    for part in dotted.split('.'):
        cur = cur[part]
    return cur


def finite(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    return x


def fmt(v) -> str:
    x = finite(v)
    if x is None:
        return 'N/A'
    if x == 0:
        return '0'
    if abs(x) < 1e-3 or abs(x) >= 1e4:
        return f'{x:.3e}'
    return f'{x:.6g}'


def relation(bilayer: float, ref: float, lower_is_better: bool = True) -> str:
    if ref == 0 and bilayer > 0:
        return 'bilayer higher than ideal zero-leakage reference'
    if bilayer == 0 and ref > 0:
        return 'zero in bilayer'
    if ref == 0 or bilayer == 0:
        return 'undefined'
    if lower_is_better:
        if bilayer < ref:
            return f'{ref / bilayer:.3g}x lower'
        return f'{bilayer / ref:.3g}x higher'
    if bilayer > ref:
        return f'{bilayer / ref:.3g}x higher'
    return f'{ref / bilayer:.3g}x lower'


def main() -> None:
    root = Path('/home/xiao/文档/KMC')
    out = root / 'comparison_al2o3_zno_vs_10nm_al2o3_envelope'
    out.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            'id': 'bilayer',
            'label': '4.5 nm Al2O3 / 6 nm ZnO bilayer KMC',
            'kind': 'KMC result',
            'path': root / 'kmc_bilayer_all/kmc_bilayer_corrected_200k_v2/transport_metrics/kmc_transport_metrics.json',
        },
        {
            'id': 'al2o3_dense_uv_ald_literature',
            'label': '10 nm high-quality UV-ALD Al2O3 literature reference',
            'kind': 'literature-calibrated',
            'path': root / 'al2o3_reference_envelope/dense_uv_ald_10nm_pet_mocon/transport_metrics/kmc_transport_metrics.json',
        },
        {
            'id': 'al2o3_pa_ald_20nm_context',
            'label': '20 nm PA-ALD Al2O3 literature context',
            'kind': 'literature-calibrated context',
            'path': root / 'al2o3_reference_envelope/pa_ald_20nm_pen_context/transport_metrics/kmc_transport_metrics.json',
        },
        {
            'id': 'al2o3_thermal_ald_literature',
            'label': '10 nm thermal ALD Al2O3 poorer-process reference',
            'kind': 'literature-calibrated',
            'path': root / 'al2o3_reference_envelope/thermal_ald_10nm_pet_mocon/transport_metrics/kmc_transport_metrics.json',
        },
        {
            'id': 'al2o3_staggered_nonthrough',
            'label': '10 nm Al2O3 staggered non-through defect idealization',
            'kind': 'connectivity-limited KMC geometry',
            'path': root / 'al2o3_reference_envelope/al2o3_10nm_staggered_nonthrough/transport_metrics/kmc_transport_metrics.json',
        },
        {
            'id': 'al2o3_through_pinhole_worst',
            'label': '10 nm Al2O3 through-pinhole worst-case KMC',
            'kind': 'KMC extreme leakage bound',
            'path': root / 'kmc_al2o3_10nm_200k/transport_metrics/kmc_transport_metrics.json',
        },
    ]
    for case in cases:
        case['metrics'] = load(case['path'])

    columns = [
        ('Ptrans', 'Ptrans'),
        ('FPT mean s', 'FPT_success_s.mean'),
        ('Lpath mean mm', 'path_length_success_mm.mean'),
        ('tau mean', 'tortuosity_success.mean'),
        ('D_eff m2/s', 'D_eff_m2_s.success_mean_FPT'),
        ('WVTR C1 g/m2/day', 'WVTR_g_m2_day.using_literature_sorbed_C1_6p77_kg_m3'),
        ('WVTR ideal vapor g/m2/day', 'WVTR_g_m2_day.using_ideal_vapor_delta_c'),
    ]

    rows = []
    for case in cases:
        row = {
            'id': case['id'],
            'label': case['label'],
            'kind': case['kind'],
            'source': str(case['path']),
        }
        for label, key in columns:
            row[label] = metric(case['metrics'], key)
        rows.append(row)

    with open(out / 'envelope_comparison_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    bilayer = cases[0]['metrics']
    bilayer_wvtr = finite(metric(bilayer, 'WVTR_g_m2_day.using_literature_sorbed_C1_6p77_kg_m3'))
    bilayer_deff = finite(metric(bilayer, 'D_eff_m2_s.success_mean_FPT'))
    reference_relations = []
    for case in cases[1:]:
        m = case['metrics']
        ref_wvtr = finite(metric(m, 'WVTR_g_m2_day.using_literature_sorbed_C1_6p77_kg_m3'))
        ref_deff = finite(metric(m, 'D_eff_m2_s.success_mean_FPT'))
        reference_relations.append({
            'reference': case['label'],
            'wvtr_relation': relation(bilayer_wvtr, ref_wvtr),
            'deff_relation': relation(bilayer_deff, ref_deff) if ref_deff is not None else 'N/A',
        })

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    plot_cases = [c for c in cases if finite(metric(c['metrics'], 'WVTR_g_m2_day.using_literature_sorbed_C1_6p77_kg_m3')) is not None]
    labels = [
        'Bilayer',
        '10 nm UV-ALD\nlit.',
        '20 nm PA-ALD\ncontext',
        '10 nm thermal\nlit.',
        'Staggered\nnon-through',
        'Through-pinhole\nworst',
    ]
    wvtr_vals = np.array([finite(metric(c['metrics'], 'WVTR_g_m2_day.using_literature_sorbed_C1_6p77_kg_m3')) for c in plot_cases], dtype=float)
    plot_vals = np.where(wvtr_vals == 0, 1e-8, wvtr_vals)
    colors = ['#d95f02', '#1b9e77', '#66a61e', '#7570b3', '#4d4d4d', '#e7298a']

    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    bars = ax.bar(range(len(plot_vals)), plot_vals, color=colors)
    ax.set_yscale('log')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=15, ha='right')
    ax.set_ylabel('WVTR, g m-2 day-1 (C1 calibration; log scale)')
    ax.set_title('Barrier comparison envelope: bilayer vs Al2O3 references')
    ax.axhline(bilayer_wvtr, color='#d95f02', lw=1.5, ls='--', alpha=0.8)
    ax.grid(axis='y', which='both', alpha=0.25)
    for bar, val, raw in zip(bars, plot_vals, wvtr_vals):
        label = '0 ideal' if raw == 0 else fmt(raw)
        ax.text(bar.get_x() + bar.get_width() / 2, val * 1.25, label, ha='center', va='bottom', fontsize=8)
    fig.savefig(out / 'wvtr_envelope_comparison.png', dpi=220)
    plt.close(fig)

    md = [
        '# Bilayer vs Al2O3 reference envelope',
        '',
        'This comparison replaces the earlier single extreme through-pinhole baseline with an envelope of Al2O3 references.',
        '',
        '| Case | Type | Ptrans | FPT mean (s) | tau mean | D_eff (m2/s) | WVTR C1 (g m-2 day-1) |',
        '|---|---|---:|---:|---:|---:|---:|',
    ]
    for row in rows:
        md.append(
            f"| {row['label']} | {row['kind']} | {fmt(row['Ptrans'])} | {fmt(row['FPT mean s'])} | "
            f"{fmt(row['tau mean'])} | {fmt(row['D_eff m2/s'])} | {fmt(row['WVTR C1 g/m2/day'])} |"
        )
    md.extend([
        '',
        '## Relations to bilayer',
        '',
        '| Reference | WVTR relation | D_eff relation |',
        '|---|---:|---:|',
    ])
    for rel in reference_relations:
        md.append(f"| {rel['reference']} | {rel['wvtr_relation']} | {rel['deff_relation']} |")
    md.extend([
        '',
        '## Interpretation',
        '',
        '- The previous 6.2e3 advantage is valid only against the through-pinhole worst-case Al2O3 model.',
        '- Against a high-quality 10 nm UV-ALD literature reference, the current bilayer KMC WVTR is about 2x lower, not thousands of times lower.',
        '- Against a poorer 10 nm thermal ALD literature reference, the bilayer is hundreds of times lower.',
        '- Against a perfect non-through-defect idealization, any finite bilayer leakage is higher than the zero-leakage reference.',
        '- Therefore the defensible claim is an envelope: bilayer performance is strongly better than through-defect or poorer-process Al2O3, comparable to or moderately better than high-quality 10 nm UV-ALD Al2O3 under the chosen calibration, but not better than an ideal fully non-percolating dense film.',
    ])
    (out / 'envelope_comparison_report.md').write_text('\n'.join(md) + '\n', encoding='utf-8')

    summary = {
        'rows': rows,
        'relations_to_bilayer': reference_relations,
        'output_dir': str(out),
    }
    (out / 'envelope_comparison_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
