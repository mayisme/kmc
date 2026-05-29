#!/usr/bin/env python3
"""Corrected Al2O3/ZnO bilayer KMC model.

Key corrections vs previous version:
1. Al2O3 dense matrix: IMPERMEABLE (no jumps allowed, not just high barrier)
2. ZnO crystalline bulk: IMPERMEABLE (wurtzite grains block water vapor)
3. ZnO grain boundaries: only diffusion pathway, but NOT straight-through
   - Introduce dead-end segments (~20% of GB blocked)
   - Closed polygonal Voronoi network forces tortuous lateral diffusion
4. Al2O3 surface: particles not hitting pinhole entrance are REFLECTED
5. D_eff formula: L^2 / (2*tau) for 1D first-passage (not 6*tau)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
from numba import cuda, int64

KB_EV_K = 8.617333262145e-5

# Cell types
CELL_AL_MATRIX = 0    # impermeable
CELL_AL_PINHOLE = 1   # fast channel
CELL_ZNO_BULK = 2     # impermeable
CELL_ZNO_GB = 3       # diffusion channel
CELL_INTERFACE = 4    # transition layer
CELL_BLOCKED = 5      # dead-end GB (impermeable)


@cuda.jit(device=True)
def xorshift64star(state):
    x = state[0]
    x ^= (x >> 12)
    x ^= (x << 25)
    x ^= (x >> 27)
    state[0] = x
    return x * 2685821657736338717


@cuda.jit(device=True)
def rand01(state):
    r = xorshift64star(state)
    val = ((r >> 11) & ((1 << 53) - 1)) * (1.0 / 9007199254740992.0)
    if val <= 1e-16:
        return 1e-16
    if val >= 1.0:
        return 0.9999999999999999
    return val


@cuda.jit(device=True)
def get_rate(cell_type, nu, inv_kbt, e_pinhole, e_gb, e_interface):
    """Return jump rate for target cell. Impermeable cells return 0."""
    if cell_type == CELL_AL_MATRIX:
        return 0.0  # impermeable
    if cell_type == CELL_ZNO_BULK:
        return 0.0  # impermeable
    if cell_type == CELL_BLOCKED:
        return 0.0  # dead-end GB
    if cell_type == CELL_AL_PINHOLE:
        return nu * math.exp(-e_pinhole * inv_kbt)
    if cell_type == CELL_ZNO_GB:
        return nu * math.exp(-e_gb * inv_kbt)
    # CELL_INTERFACE
    return nu * math.exp(-e_interface * inv_kbt)


@cuda.jit
def bilayer_kernel(grid, width, height, max_events, max_time_s, nu, inv_kbt,
                   e_pinhole, e_gb, e_interface, e_surface,
                   pinhole_mask, al_thickness,
                   seeds, first_times, steps, transmitted, final_y):
    """Each thread simulates one particle.
    
    Surface diffusion is handled analytically: particle starts at random x,
    does 1D random walk on periodic surface until hitting a pinhole.
    The expected time is computed from the surface hop rate and distance.
    Then bulk transport through the structure is simulated with full kMC.
    """
    i = cuda.grid(1)
    if i >= first_times.size:
        return

    st = cuda.local.array(1, dtype=int64)
    st[0] = seeds[i]
    x = int(rand01(st) * width)
    if x >= width:
        x = width - 1

    # --- Surface diffusion phase (simulated step by step on 1D surface) ---
    # Find nearest pinhole by random walk on surface
    k_surf = nu * math.exp(-e_surface * inv_kbt)
    t_acc = 0.0
    surf_steps = 0
    max_surf_steps = max_events // 2  # reserve half events for bulk

    while pinhole_mask[x] == 0 and surf_steps < max_surf_steps:
        # Two choices: left or right (equal rate)
        dt = -math.log(rand01(st)) / (2.0 * k_surf)
        t_acc += dt
        if t_acc > max_time_s:
            steps[i] = surf_steps
            first_times[i] = t_acc
            transmitted[i] = 0
            final_y[i] = 0
            seeds[i] = st[0]
            return
        if rand01(st) < 0.5:
            x = (x - 1) % width
        else:
            x = (x + 1) % width
        surf_steps += 1

    if pinhole_mask[x] == 0:
        # Did not find pinhole in allowed steps
        steps[i] = surf_steps
        first_times[i] = t_acc
        transmitted[i] = 0
        final_y[i] = 0
        seeds[i] = st[0]
        return

    # --- Bulk transport phase: start at y=0 (pinhole entrance) ---
    y = 0
    bulk_events = max_events - surf_steps

    for ev in range(bulk_events):
        xl = (x - 1) % width
        xr = (x + 1) % width
        yu = y - 1
        yd = y + 1

        # Up
        k_up = 0.0
        if yu >= 0:
            c_up = grid[yu, x]
            k_up = get_rate(c_up, nu, inv_kbt, e_pinhole, e_gb, e_interface)
        # y=0 up goes to surface - allow return with surface rate
        # (but we don't simulate surface again, just bounce back)

        # Down
        k_down = 0.0
        if yd < height:
            c_down = grid[yd, x]
            k_down = get_rate(c_down, nu, inv_kbt, e_pinhole, e_gb, e_interface)
        elif yd >= height:
            k_down = nu * math.exp(-e_gb * inv_kbt)

        # Left
        c_left = grid[y, xl]
        k_left = get_rate(c_left, nu, inv_kbt, e_pinhole, e_gb, e_interface)

        # Right
        c_right = grid[y, xr]
        k_right = get_rate(c_right, nu, inv_kbt, e_pinhole, e_gb, e_interface)

        k_sum = k_up + k_down + k_left + k_right
        if k_sum <= 0.0:
            steps[i] = surf_steps + ev
            first_times[i] = t_acc
            transmitted[i] = 0
            final_y[i] = y
            seeds[i] = st[0]
            return

        t_acc += -math.log(rand01(st)) / k_sum
        if t_acc > max_time_s:
            steps[i] = surf_steps + ev + 1
            first_times[i] = t_acc
            transmitted[i] = 0
            final_y[i] = y
            seeds[i] = st[0]
            return

        u = rand01(st) * k_sum
        if u < k_up:
            if yu >= 0:
                y = yu
            # else: bounce (stay at y=0)
        elif u < k_up + k_down:
            if yd >= height:
                steps[i] = surf_steps + ev + 1
                first_times[i] = t_acc
                transmitted[i] = 1
                final_y[i] = height
                seeds[i] = st[0]
                return
            y = yd
        elif u < k_up + k_down + k_left:
            x = xl
        else:
            x = xr

    steps[i] = max_events
    first_times[i] = t_acc
    transmitted[i] = 0
    final_y[i] = y
    seeds[i] = st[0]



def draw_disk(grid, y0, x0, radius_cells, value):
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


def build_bilayer_structure(width_nm=320.0, dx_nm=0.5, seed=20260515,
                            periods=1, al_nm=4.5, zno_nm=6.0,
                            pinhole_spacing_min_nm=50.0, pinhole_spacing_max_nm=80.0,
                            pinhole_d_min_nm=1.0, pinhole_d_max_nm=2.0,
                            grain_min_nm=5.0, grain_max_nm=10.0,
                            offset_min_nm=25.0, offset_max_nm=50.0,
                            gb_block_fraction=0.20):
    """Build corrected bilayer structure.
    
    Key differences:
    - Al2O3 matrix cells are CELL_AL_MATRIX (impermeable, rate=0)
    - ZnO bulk cells are CELL_ZNO_BULK (impermeable, rate=0)
    - ZnO GB network from Voronoi, with gb_block_fraction of GB cells
      converted to CELL_BLOCKED (dead ends, rate=0)
    - Returns pinhole_mask for surface reflection
    """
    rng = np.random.default_rng(seed)
    width = int(round(width_nm / dx_nm))
    al_cells = int(round(al_nm / dx_nm))
    zno_cells = int(round(zno_nm / dx_nm))
    period_cells = al_cells + zno_cells
    height = periods * period_cells
    grid = np.full((height, width), CELL_AL_MATRIX, dtype=np.int32)
    # Pinhole entrance mask at y=0 (1=pinhole, 0=reflected)
    pinhole_mask = np.zeros(width, dtype=np.int32)
    metadata = {'periods': []}

    prev_centers = None
    for p in range(periods):
        y0 = p * period_cells
        al_start, al_end = y0, y0 + al_cells
        z_start, z_end = al_end, y0 + period_cells

        # Fill ZnO region with bulk (impermeable)
        grid[z_start:z_end, :] = CELL_ZNO_BULK
        # Interface row: BOTH the last Al row AND first ZnO row are INTERFACE
        # This allows particles exiting pinhole to diffuse laterally to nearest GB
        grid[al_end - 1, :] = CELL_INTERFACE
        grid[z_start, :] = CELL_INTERFACE

        # --- Al2O3 pinholes ---
        centers = []
        if prev_centers is None:
            x_nm = rng.uniform(0, pinhole_spacing_max_nm)
            while x_nm < width_nm:
                centers.append(x_nm)
                x_nm += rng.uniform(pinhole_spacing_min_nm, pinhole_spacing_max_nm)
        else:
            sign = -1 if p % 2 else 1
            offset = sign * rng.uniform(offset_min_nm, offset_max_nm)
            centers = sorted([((c + offset) % width_nm) for c in prev_centers])
        prev_centers = centers

        pinholes = []
        for c_nm in centers:
            d_nm = rng.uniform(pinhole_d_min_nm, pinhole_d_max_nm)
            r = max(1.0, (d_nm / 2.0) / dx_nm)
            cx = c_nm / dx_nm
            for yy in range(al_start, al_end):
                jitter = rng.normal(0.0, 0.25)
                draw_disk(grid, yy, cx + jitter, r, CELL_AL_PINHOLE)
            # Mark pinhole entrance in mask (at y=0 for first period)
            if p == 0:
                cx_int = int(round(cx)) % width
                r_int = max(1, int(round(r)))
                for dx_i in range(-r_int, r_int + 1):
                    xi = (cx_int + dx_i) % width
                    pinhole_mask[xi] = 1
            pinholes.append({'x_nm': float(c_nm), 'diameter_nm': float(d_nm)})

        # --- ZnO Voronoi grain boundary network ---
        zno_h = z_end - z_start
        mean_grain_nm = 0.5 * (grain_min_nm + grain_max_nm)
        mean_grain_cells = max(2.0, mean_grain_nm / dx_nm)
        zno_area_cells = width * zno_h
        n_grains = max(4, int(round(zno_area_cells / (mean_grain_cells ** 2))))

        nuclei_x = rng.uniform(0, width, size=n_grains)
        nuclei_y = rng.uniform(0, zno_h, size=n_grains)

        # Offset nuclei from pinholes to prevent vertical alignment
        for c_nm in centers:
            off_nm = rng.choice([-1.0, 1.0]) * rng.uniform(offset_min_nm, offset_max_nm)
            nuclei_x = np.append(nuclei_x, ((c_nm + off_nm) / dx_nm) % width)
            nuclei_y = np.append(nuclei_y, rng.uniform(0, zno_h))
        n_grains = nuclei_x.size

        # Voronoi labeling
        labels = np.empty((zno_h, width), dtype=np.int32)
        for yy in range(zno_h):
            for xx in range(width):
                best = 0
                best_d2 = 1.0e30
                for gi in range(n_grains):
                    dxp = abs(xx - nuclei_x[gi])
                    if dxp > 0.5 * width:
                        dxp = width - dxp
                    dyp = yy - nuclei_y[gi]
                    # Slight columnar anisotropy
                    d2 = dxp * dxp + 0.65 * dyp * dyp
                    if d2 < best_d2:
                        best_d2 = d2
                        best = gi
                labels[yy, xx] = best

        # Extract GB mask
        gb_mask = np.zeros((zno_h, width), dtype=bool)
        for yy in range(zno_h):
            for xx in range(width):
                lab = labels[yy, xx]
                if labels[yy, (xx - 1) % width] != lab or labels[yy, (xx + 1) % width] != lab:
                    gb_mask[yy, xx] = True
                elif yy > 0 and labels[yy - 1, xx] != lab:
                    gb_mask[yy, xx] = True
                elif yy + 1 < zno_h and labels[yy + 1, xx] != lab:
                    gb_mask[yy, xx] = True

        # Mark GB cells
        grid[z_start:z_end, :][gb_mask] = CELL_ZNO_GB

        # Restore interface rows (first and last ZnO row) after GB marking
        # First ZnO row = lateral diffusion layer from pinhole exit to GB
        grid[z_start, :] = CELL_INTERFACE
        # Last ZnO row = exit interface
        grid[z_end - 1, :] = CELL_INTERFACE

        # --- Introduce dead-end blockages in GB network ---
        # Randomly block gb_block_fraction of GB cells to prevent straight-through
        # But preserve first and last rows for entry/exit connectivity
        gb_coords = np.argwhere(gb_mask)
        # Only block interior GB cells (not row 0 or row zno_h-1)
        interior_mask = (gb_coords[:, 0] > 0) & (gb_coords[:, 0] < zno_h - 1)
        interior_coords = gb_coords[interior_mask]
        n_block = int(round(len(interior_coords) * gb_block_fraction))
        if n_block > 0:
            block_idx = rng.choice(len(interior_coords), size=n_block, replace=False)
            for idx in block_idx:
                by, bx = interior_coords[idx]
                grid[z_start + by, bx] = CELL_BLOCKED

        gb_frac = float(gb_mask.mean())
        metadata['periods'].append({
            'period': p + 1, 'pinholes': pinholes,
            'zno_grain_count': int(n_grains),
            'zno_gb_fraction': gb_frac,
            'gb_blocked_fraction': float(gb_block_fraction),
            'pinhole_coverage': float(pinhole_mask.sum() / width),
        })

    return grid, width, height, al_cells, zno_cells, pinhole_mask, metadata



def summarize(first_times, steps, transmitted, final_y, total_thickness_nm, pinhole_mask, width):
    ok = transmitted.astype(bool)
    n = int(transmitted.size)
    n_ok = int(ok.sum())
    # With surface diffusion, all particles enter the simulation (no reflection)
    out = {
        'particles': n,
        'transmitted': n_ok,
        'success_rate': float(n_ok / n) if n > 0 else 0.0,
        'pinhole_area_fraction': float(pinhole_mask.sum() / width),
    }
    if n_ok > 0:
        ft = first_times[ok]
        st = steps[ok]
        out.update({
            'mean_first_passage_s': float(np.mean(ft)),
            'median_first_passage_s': float(np.median(ft)),
            'p10_first_passage_s': float(np.percentile(ft, 10)),
            'p90_first_passage_s': float(np.percentile(ft, 90)),
            'mean_path_steps_success': float(np.mean(st)),
            'median_path_steps_success': float(np.median(st)),
            # Corrected: 1D first-passage D = L^2 / (2*tau)
            'D_eff_m2_s': float(((total_thickness_nm * 1e-9) ** 2) / (2.0 * np.mean(ft))),
            'tortuosity_steps_over_thickness': float(np.mean(st) / max(total_thickness_nm / 0.5, 1.0)),
        })
    else:
        for k in ['mean_first_passage_s', 'median_first_passage_s', 'p10_first_passage_s',
                   'p90_first_passage_s', 'mean_path_steps_success', 'median_path_steps_success',
                   'D_eff_m2_s', 'tortuosity_steps_over_thickness']:
            out[k] = None
    # Non-transmitted analysis
    not_trans = ~ok
    if not_trans.sum() > 0:
        fy = final_y[not_trans] * 0.5  # convert to nm
        out['non_transmitted_trapped'] = int(not_trans.sum())
        out['trapped_mean_depth_nm'] = float(np.mean(fy))
        out['trapped_median_depth_nm'] = float(np.median(fy))
    return out


def save_summary_png(path, grid, result, args, pinhole_mask):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    # Structure map
    cmap = ListedColormap(['#08306b', '#41b6c4', '#fed976', '#f03b20', '#7a0177', '#333333'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    im = axes[0, 0].imshow(grid, aspect='auto', origin='upper', cmap=cmap, norm=norm)
    axes[0, 0].set_title(f'Corrected bilayer structure (periods={args.periods})')
    axes[0, 0].set_xlabel('x cell')
    axes[0, 0].set_ylabel('z cell')
    cbar = fig.colorbar(im, ax=axes[0, 0], ticks=[0, 1, 2, 3, 4, 5])
    cbar.ax.set_yticklabels(['Al matrix\n(imper.)', 'Al pinhole', 'ZnO bulk\n(imper.)',
                              'ZnO GB', 'interface', 'blocked GB'])

    # Barrier info
    labels = ['Al matrix', 'Al pinhole', 'ZnO bulk', 'ZnO GB', 'interface', 'blocked GB']
    vals = ['∞ (imper.)', f'{args.e_al_pinhole} eV', '∞ (imper.)',
            f'{args.e_zno_gb} eV', f'{args.e_interface} eV', '∞ (imper.)']
    axes[0, 1].barh(labels, [0, args.e_al_pinhole, 0, args.e_zno_gb, args.e_interface, 0],
                    color=['#08306b', '#41b6c4', '#fed976', '#f03b20', '#7a0177', '#333333'])
    axes[0, 1].set_xlabel('E_ij (eV), 0 = impermeable')
    axes[0, 1].set_title('Migration barriers (corrected)')

    # FPT distribution
    data = np.load(path / 'bilayer_corrected_raw.npz')
    ft = data['first_times']
    trans = data['transmitted'].astype(bool)
    if trans.sum() > 0:
        axes[1, 0].hist(ft[trans], bins=60, color='#3182bd', alpha=0.85)
        axes[1, 0].set_xlabel('First-passage time (s)')
        axes[1, 0].set_ylabel('count')
        axes[1, 0].set_title('Successful particle FPT distribution')
    else:
        axes[1, 0].text(0.5, 0.5, 'No successful transmission', ha='center', va='center',
                        fontsize=14, color='red')
        axes[1, 0].set_axis_off()

    # Summary text
    def fmt(v):
        if v is None:
            return 'N/A'
        if isinstance(v, float) and (abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0)):
            return f'{v:.4g}'
        if isinstance(v, float):
            return f'{v:.4f}'
        return str(v)

    txt = '\n'.join([
        f"particles = {result['particles']}",
        f"transmitted = {result['transmitted']}",
        f"success_rate = {fmt(result['success_rate'])}",
        f"pinhole area frac = {fmt(result['pinhole_area_fraction'])}",
        f"mean FPT = {fmt(result.get('mean_first_passage_s'))} s",
        f"median FPT = {fmt(result.get('median_first_passage_s'))} s",
        f"D_eff = {fmt(result.get('D_eff_m2_s'))} m²/s",
        f"elapsed = {result['elapsed_s']:.2f} s",
    ])
    axes[1, 1].text(0.03, 0.97, txt, va='top', ha='left', family='monospace', fontsize=11)
    axes[1, 1].set_axis_off()
    axes[1, 1].set_title('Run summary (corrected model)')
    fig.savefig(path / 'bilayer_corrected_summary.png', dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Corrected Al2O3/ZnO bilayer KMC')
    ap.add_argument('--particles', type=int, default=200000)
    ap.add_argument('--max-events', type=int, default=200000)
    ap.add_argument('--max-time-s', type=float, default=1e10)
    ap.add_argument('--width-nm', type=float, default=320.0)
    ap.add_argument('--dx-nm', type=float, default=0.5)
    ap.add_argument('--periods', type=int, default=1)
    ap.add_argument('--threads-per-block', type=int, default=256)
    ap.add_argument('--out-dir', type=str, default='kmc_bilayer_corrected')
    ap.add_argument('--seed', type=int, default=20260515)
    ap.add_argument('--temperature-k', type=float, default=311.15)
    ap.add_argument('--rh', type=float, default=0.90)
    ap.add_argument('--nu', type=float, default=1.0e12)
    # Only permeable cells need barriers
    ap.add_argument('--e-al-pinhole', type=float, default=0.58)
    ap.add_argument('--e-zno-gb', type=float, default=0.72,
                    help='GB barrier raised from 0.64 to 0.72 eV for more realistic tortuous diffusion')
    ap.add_argument('--e-interface', type=float, default=0.62,
                    help='Interface barrier: between pinhole and GB, allows lateral diffusion to find GB entry')
    ap.add_argument('--e-surface', type=float, default=0.45,
                    help='Surface adsorption diffusion barrier on Al2O3 (90%RH adsorbed water layer)')
    ap.add_argument('--gb-block-fraction', type=float, default=0.20,
                    help='Fraction of GB cells blocked to prevent straight-through paths')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not cuda.is_available():
        raise SystemExit('Numba CUDA is not available')

    t0 = time.time()
    grid, width, height, al_cells, zno_cells, pinhole_mask, geom_meta = build_bilayer_structure(
        args.width_nm, args.dx_nm, args.seed, periods=args.periods,
        gb_block_fraction=args.gb_block_fraction)

    rng = np.random.default_rng(args.seed + 777)
    seeds = rng.integers(1, np.iinfo(np.int64).max, size=args.particles, dtype=np.int64)
    first_times = np.zeros(args.particles, dtype=np.float64)
    steps = np.zeros(args.particles, dtype=np.int32)
    transmitted = np.zeros(args.particles, dtype=np.int32)
    final_y = np.zeros(args.particles, dtype=np.int32)

    d_grid = cuda.to_device(grid)
    d_pinhole_mask = cuda.to_device(pinhole_mask)
    d_seeds = cuda.to_device(seeds)
    d_first = cuda.to_device(first_times)
    d_steps = cuda.to_device(steps)
    d_trans = cuda.to_device(transmitted)
    d_final_y = cuda.to_device(final_y)

    inv_kbt = 1.0 / (KB_EV_K * args.temperature_k)
    blocks = (args.particles + args.threads_per_block - 1) // args.threads_per_block
    bilayer_kernel[blocks, args.threads_per_block](
        d_grid, width, height, args.max_events, args.max_time_s, args.nu, inv_kbt,
        args.e_al_pinhole, args.e_zno_gb, args.e_interface, args.e_surface,
        d_pinhole_mask, al_cells,
        d_seeds, d_first, d_steps, d_trans, d_final_y)
    cuda.synchronize()

    first_times = d_first.copy_to_host()
    steps = d_steps.copy_to_host()
    transmitted = d_trans.copy_to_host()
    final_y = d_final_y.copy_to_host()
    elapsed = time.time() - t0

    total_thickness_nm = args.periods * (args.al_nm if hasattr(args, 'al_nm') else 4.5 + 6.0)
    total_thickness_nm = args.periods * (4.5 + 6.0)
    result = summarize(first_times, steps, transmitted, final_y, total_thickness_nm, pinhole_mask, width)
    result.update({
        'elapsed_s': float(elapsed),
        'model': 'corrected_bilayer_v2',
        'corrections': [
            'Al2O3_matrix_impermeable',
            'ZnO_bulk_impermeable',
            'surface_reflection_mechanism',
            'GB_dead_end_blockage',
            'D_eff_1D_formula',
        ],
        'environment': {'temperature_c': 38.0, 'temperature_k': args.temperature_k, 'rh': args.rh},
        'periods': int(args.periods),
        'total_thickness_nm': total_thickness_nm,
        'width_nm': args.width_nm,
        'dx_nm': args.dx_nm,
        'max_events': int(args.max_events),
        'barriers_eV': {
            'E_Al_pinhole': args.e_al_pinhole,
            'E_ZnO_GB': args.e_zno_gb,
            'E_interface': args.e_interface,
            'E_surface': args.e_surface,
            'E_Al_matrix': 'infinite (impermeable)',
            'E_ZnO_bulk': 'infinite (impermeable)',
        },
        'gb_block_fraction': args.gb_block_fraction,
        'geometry_metadata': geom_meta,
    })

    with open(out_dir / 'bilayer_corrected_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['metric', 'value'])
        def emit(prefix, obj):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    emit(f'{prefix}.{k}' if prefix else k, v)
            elif isinstance(obj, list):
                w.writerow([prefix, json.dumps(obj, ensure_ascii=False)])
            else:
                w.writerow([prefix, obj])
        for k, v in result.items():
            emit(k, v)

    with open(out_dir / 'bilayer_corrected_results.json', 'w') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    np.savez_compressed(out_dir / 'bilayer_corrected_raw.npz',
                        first_times=first_times, steps=steps, transmitted=transmitted,
                        final_y=final_y, grid=grid, pinhole_mask=pinhole_mask)
    save_summary_png(out_dir, grid, result, args, pinhole_mask)
    print(json.dumps({k: v for k, v in result.items() if k != 'geometry_metadata'}, indent=2, ensure_ascii=False))
    print(f'outputs: {out_dir.resolve()}')


if __name__ == '__main__':
    main()
