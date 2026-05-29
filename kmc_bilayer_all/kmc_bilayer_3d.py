#!/usr/bin/env python3
"""3D on-lattice Al2O3/ZnO bilayer kMC model.

Geometry convention: grid[z, x, y]
- z is physical thickness: 4.5 nm Al2O3 + 6.0 nm ZnO by default.
- x/y are lateral periodic dimensions.
- Al2O3 and ZnO bulk are impermeable; Al pinholes, interfaces, and ZnO GBs are permeable.

This script intentionally does not replace kmc_bilayer_corrected.py. It is a 3D extension
for thin-film lateral-mismatch statistics and tortuous paths.
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

CELL_AL_MATRIX = 0
CELL_AL_PINHOLE = 1
CELL_ZNO_BULK = 2
CELL_ZNO_GB = 3
CELL_INTERFACE = 4
CELL_BLOCKED = 5


@cuda.jit(device=True)
def xorshift64star(state):
    x = state[0]
    x ^= x >> 12
    x ^= x << 25
    x ^= x >> 27
    state[0] = x
    return x * 2685821657736338717


@cuda.jit(device=True)
def rand01(state):
    r = xorshift64star(state)
    val = ((r >> 11) & ((1 << 53) - 1)) * (1.0 / 9007199254740992.0)
    if val <= 1.0e-16:
        return 1.0e-16
    if val >= 1.0:
        return 0.9999999999999999
    return val


@cuda.jit(device=True)
def get_rate(cell_type, nu, inv_kbt, e_pinhole, e_gb, e_interface):
    if cell_type == CELL_AL_MATRIX:
        return 0.0
    if cell_type == CELL_ZNO_BULK:
        return 0.0
    if cell_type == CELL_BLOCKED:
        return 0.0
    if cell_type == CELL_AL_PINHOLE:
        return nu * math.exp(-e_pinhole * inv_kbt)
    if cell_type == CELL_ZNO_GB:
        return nu * math.exp(-e_gb * inv_kbt)
    return nu * math.exp(-e_interface * inv_kbt)


@cuda.jit
def bilayer_3d_kernel(grid, nx, ny, nz, max_events, max_time_s, nu, inv_kbt,
                      e_pinhole, e_gb, e_interface, e_surface, pinhole_mask,
                      seeds, first_times, steps, transmitted, final_z, final_x, final_y):
    i = cuda.grid(1)
    if i >= first_times.size:
        return

    st = cuda.local.array(1, dtype=int64)
    st[0] = seeds[i]

    x = int(rand01(st) * nx)
    y = int(rand01(st) * ny)
    if x >= nx:
        x = nx - 1
    if y >= ny:
        y = ny - 1

    k_surf = nu * math.exp(-e_surface * inv_kbt)
    t_acc = 0.0
    surf_steps = 0
    max_surf_steps = max_events // 2

    # 2D periodic surface diffusion until a pinhole mouth is found.
    while pinhole_mask[x, y] == 0 and surf_steps < max_surf_steps:
        t_acc += -math.log(rand01(st)) / (4.0 * k_surf)
        if t_acc > max_time_s:
            steps[i] = surf_steps
            first_times[i] = t_acc
            transmitted[i] = 0
            final_z[i] = 0
            final_x[i] = x
            final_y[i] = y
            seeds[i] = st[0]
            return
        r = rand01(st)
        if r < 0.25:
            x = (x - 1) % nx
        elif r < 0.50:
            x = (x + 1) % nx
        elif r < 0.75:
            y = (y - 1) % ny
        else:
            y = (y + 1) % ny
        surf_steps += 1

    if pinhole_mask[x, y] == 0:
        steps[i] = surf_steps
        first_times[i] = t_acc
        transmitted[i] = 0
        final_z[i] = 0
        final_x[i] = x
        final_y[i] = y
        seeds[i] = st[0]
        return

    z = 0
    bulk_events = max_events - surf_steps
    for ev in range(bulk_events):
        zm = z - 1
        zp = z + 1
        xm = (x - 1) % nx
        xp = (x + 1) % nx
        ym = (y - 1) % ny
        yp = (y + 1) % ny

        k_zm = 0.0
        if zm >= 0:
            k_zm = get_rate(grid[zm, x, y], nu, inv_kbt, e_pinhole, e_gb, e_interface)
        else:
            k_zm = k_surf

        k_zp = 0.0
        if zp < nz:
            k_zp = get_rate(grid[zp, x, y], nu, inv_kbt, e_pinhole, e_gb, e_interface)
        else:
            k_zp = nu * math.exp(-e_gb * inv_kbt)

        k_xm = get_rate(grid[z, xm, y], nu, inv_kbt, e_pinhole, e_gb, e_interface)
        k_xp = get_rate(grid[z, xp, y], nu, inv_kbt, e_pinhole, e_gb, e_interface)
        k_ym = get_rate(grid[z, x, ym], nu, inv_kbt, e_pinhole, e_gb, e_interface)
        k_yp = get_rate(grid[z, x, yp], nu, inv_kbt, e_pinhole, e_gb, e_interface)

        k_sum = k_zm + k_zp + k_xm + k_xp + k_ym + k_yp
        if k_sum <= 0.0:
            steps[i] = surf_steps + ev
            first_times[i] = t_acc
            transmitted[i] = 0
            final_z[i] = z
            final_x[i] = x
            final_y[i] = y
            seeds[i] = st[0]
            return

        t_acc += -math.log(rand01(st)) / k_sum
        if t_acc > max_time_s:
            steps[i] = surf_steps + ev + 1
            first_times[i] = t_acc
            transmitted[i] = 0
            final_z[i] = z
            final_x[i] = x
            final_y[i] = y
            seeds[i] = st[0]
            return

        u = rand01(st) * k_sum
        if u < k_zm:
            if zm >= 0:
                z = zm
        elif u < k_zm + k_zp:
            if zp >= nz:
                steps[i] = surf_steps + ev + 1
                first_times[i] = t_acc
                transmitted[i] = 1
                final_z[i] = nz
                final_x[i] = x
                final_y[i] = y
                seeds[i] = st[0]
                return
            z = zp
        elif u < k_zm + k_zp + k_xm:
            x = xm
        elif u < k_zm + k_zp + k_xm + k_xp:
            x = xp
        elif u < k_zm + k_zp + k_xm + k_xp + k_ym:
            y = ym
        else:
            y = yp

    steps[i] = max_events
    first_times[i] = t_acc
    transmitted[i] = 0
    final_z[i] = z
    final_x[i] = x
    final_y[i] = y
    seeds[i] = st[0]


def _periodic_delta(a, b, n):
    d = abs(a - b)
    return min(d, n - d)


def draw_cylinder_xy(grid, z0, z1, cx, cy, radius_cells, value):
    nz, nx, ny = grid.shape
    r2 = radius_cells * radius_cells
    xmin = int(math.floor(cx - radius_cells - 1))
    xmax = int(math.ceil(cx + radius_cells + 2))
    ymin = int(math.floor(cy - radius_cells - 1))
    ymax = int(math.ceil(cy + radius_cells + 2))
    for z in range(max(0, z0), min(nz, z1)):
        for xx in range(xmin, xmax):
            x = xx % nx
            dx = _periodic_delta(x, cx % nx, nx)
            for yy in range(ymin, ymax):
                y = yy % ny
                dy = _periodic_delta(y, cy % ny, ny)
                if dx * dx + dy * dy <= r2:
                    grid[z, x, y] = value


def draw_sphere(grid, cz, cx, cy, radius_cells, value):
    nz, nx, ny = grid.shape
    r2 = radius_cells * radius_cells
    for z in range(max(0, int(math.floor(cz - radius_cells - 1))), min(nz, int(math.ceil(cz + radius_cells + 2)))):
        dz = z - cz
        for xx in range(int(math.floor(cx - radius_cells - 1)), int(math.ceil(cx + radius_cells + 2))):
            x = xx % nx
            dx = _periodic_delta(x, cx % nx, nx)
            for yy in range(int(math.floor(cy - radius_cells - 1)), int(math.ceil(cy + radius_cells + 2))):
                y = yy % ny
                dy = _periodic_delta(y, cy % ny, ny)
                if dx * dx + dy * dy + dz * dz <= r2:
                    grid[z, x, y] = value


def build_bilayer_3d_structure(width_x_nm=120.0, width_y_nm=120.0, dx_nm=0.5, seed=20260529,
                               al_nm=4.5, zno_nm=6.0,
                               pinhole_spacing_min_nm=50.0, pinhole_spacing_max_nm=80.0,
                               pinhole_d_min_nm=1.0, pinhole_d_max_nm=2.0,
                               al_free_volume_count=18,
                               al_free_volume_d_min_nm=0.8, al_free_volume_d_max_nm=1.8,
                               grain_min_nm=5.0, grain_max_nm=10.0,
                               zno_gb_width_cells=1,
                               zno_gb_wander_nm=1.0,
                               offset_min_nm=20.0, offset_max_nm=50.0,
                               gb_block_fraction=0.20):
    rng = np.random.default_rng(seed)
    nx = int(round(width_x_nm / dx_nm))
    ny = int(round(width_y_nm / dx_nm))
    al_cells = int(round(al_nm / dx_nm))
    zno_cells = int(round(zno_nm / dx_nm))
    nz = al_cells + zno_cells

    grid = np.full((nz, nx, ny), CELL_AL_MATRIX, dtype=np.int32)
    grid[al_cells:nz, :, :] = CELL_ZNO_BULK
    grid[al_cells - 1, :, :] = CELL_INTERFACE
    grid[al_cells, :, :] = CELL_INTERFACE
    grid[nz - 1, :, :] = CELL_INTERFACE
    pinhole_mask = np.zeros((nx, ny), dtype=np.int32)

    pinholes = []
    mean_spacing = 0.5 * (pinhole_spacing_min_nm + pinhole_spacing_max_nm)
    area_nm2 = width_x_nm * width_y_nm
    n_pinholes = max(1, int(round(area_nm2 / (mean_spacing * mean_spacing))))
    for _ in range(n_pinholes):
        cx_nm = rng.uniform(0.0, width_x_nm)
        cy_nm = rng.uniform(0.0, width_y_nm)
        d_nm = rng.uniform(pinhole_d_min_nm, pinhole_d_max_nm)
        r = max(1.0, (d_nm / 2.0) / dx_nm)
        cx = cx_nm / dx_nm
        cy = cy_nm / dx_nm
        for z in range(al_cells):
            drift_x = rng.normal(0.0, 0.20)
            drift_y = rng.normal(0.0, 0.20)
            draw_cylinder_xy(grid, z, z + 1, cx + drift_x, cy + drift_y, r, CELL_AL_PINHOLE)
        ci = int(round(cx)) % nx
        cj = int(round(cy)) % ny
        ri = max(1, int(round(r)))
        for xx in range(ci - ri, ci + ri + 1):
            x = xx % nx
            for yy in range(cj - ri, cj + ri + 1):
                y = yy % ny
                if _periodic_delta(x, ci, nx) ** 2 + _periodic_delta(y, cj, ny) ** 2 <= ri * ri:
                    pinhole_mask[x, y] = 1
                    # Keep the surface mask and actual z=0 lattice mouth exactly consistent.
                    grid[0, x, y] = CELL_AL_PINHOLE
        pinholes.append({'x_nm': float(cx_nm), 'y_nm': float(cy_nm), 'diameter_nm': float(d_nm)})

    free_volume_defects = []
    for _ in range(al_free_volume_count):
        d_nm = rng.uniform(al_free_volume_d_min_nm, al_free_volume_d_max_nm)
        r = max(1.0, (d_nm / 2.0) / dx_nm)
        if al_cells > 4:
            cz = rng.uniform(1.0, al_cells - 2.0)
        else:
            cz = max(0.0, 0.5 * (al_cells - 1))
        cx = rng.uniform(0, nx)
        cy = rng.uniform(0, ny)
        draw_sphere(grid, cz, cx, cy, r, CELL_AL_PINHOLE)
        # Do not allow local free-volume pockets to become surface mouths.
        grid[0][pinhole_mask == 0] = np.where(grid[0][pinhole_mask == 0] == CELL_AL_PINHOLE,
                                               CELL_AL_MATRIX, grid[0][pinhole_mask == 0])
        grid[al_cells - 1][pinhole_mask == 0] = np.where(grid[al_cells - 1][pinhole_mask == 0] == CELL_AL_PINHOLE,
                                                          CELL_INTERFACE, grid[al_cells - 1][pinhole_mask == 0])
        free_volume_defects.append({'x_nm': float(cx * dx_nm), 'y_nm': float(cy * dx_nm),
                                    'diameter_nm': float(d_nm), 'through_pinhole': False})

    zno_h = zno_cells
    mean_grain_nm = 0.5 * (grain_min_nm + grain_max_nm)
    cell_area = max(1.0, (mean_grain_nm / dx_nm) ** 2)
    n_grains = max(4, int(round((nx * ny) / cell_area)))
    nuclei_x = rng.uniform(0, nx, size=n_grains)
    nuclei_y = rng.uniform(0, ny, size=n_grains)

    # Add GB-attracting nuclei laterally offset from pinholes so pinholes seldom align directly with GB mouths.
    for p in pinholes:
        theta = rng.uniform(0, 2.0 * math.pi)
        off = rng.uniform(offset_min_nm, offset_max_nm) / dx_nm
        nuclei_x = np.append(nuclei_x, (p['x_nm'] / dx_nm + off * math.cos(theta)) % nx)
        nuclei_y = np.append(nuclei_y, (p['y_nm'] / dx_nm + off * math.sin(theta)) % ny)
    n_grains = int(nuclei_x.size)

    labels = np.empty((nx, ny), dtype=np.int32)
    for x in range(nx):
        for y in range(ny):
            best = 0
            best_d2 = 1.0e30
            for gi in range(n_grains):
                dxp = _periodic_delta(x, nuclei_x[gi], nx)
                dyp = _periodic_delta(y, nuclei_y[gi], ny)
                d2 = dxp * dxp + dyp * dyp
                if d2 < best_d2:
                    best_d2 = d2
                    best = gi
            labels[x, y] = best

    gb2d = np.zeros((nx, ny), dtype=bool)
    for x in range(nx):
        for y in range(ny):
            lab = labels[x, y]
            if labels[(x - 1) % nx, y] != lab or labels[(x + 1) % nx, y] != lab:
                gb2d[x, y] = True
            elif labels[x, (y - 1) % ny] != lab or labels[x, (y + 1) % ny] != lab:
                gb2d[x, y] = True

    gb3d = np.zeros((zno_h, nx, ny), dtype=bool)
    amp = max(0.0, zno_gb_wander_nm / dx_nm)
    for zz in range(zno_h):
        shift_x = int(round(amp * math.sin(2.0 * math.pi * zz / max(1, zno_h - 1))))
        shift_y = int(round(amp * math.cos(2.0 * math.pi * zz / max(1, zno_h - 1))))
        shifted = np.roll(np.roll(gb2d, shift_x, axis=0), shift_y, axis=1)
        gb3d[zz] = shifted
        for w in range(1, max(1, zno_gb_width_cells)):
            gb3d[zz] |= np.roll(shifted, w, axis=0) | np.roll(shifted, -w, axis=0)
            gb3d[zz] |= np.roll(shifted, w, axis=1) | np.roll(shifted, -w, axis=1)

    z_start = al_cells
    grid[z_start:nz, :, :][gb3d] = CELL_ZNO_GB
    grid[z_start, :, :] = CELL_INTERFACE
    grid[nz - 1, :, :] = CELL_INTERFACE

    gb_coords = np.argwhere(gb3d)
    interior = gb_coords[(gb_coords[:, 0] > 0) & (gb_coords[:, 0] < zno_h - 1)]
    n_block = int(round(len(interior) * gb_block_fraction))
    if n_block > 0:
        idx = rng.choice(len(interior), size=n_block, replace=False)
        for row in interior[idx]:
            zz, x, y = int(row[0]), int(row[1]), int(row[2])
            grid[z_start + zz, x, y] = CELL_BLOCKED

    metadata = {
        'geometry': '3d_on_lattice_bilayer',
        'shape_zxy': [int(nz), int(nx), int(ny)],
        'al_cells': int(al_cells),
        'zno_cells': int(zno_cells),
        'total_thickness_nm': float(al_nm + zno_nm),
        'pinhole_count': int(len(pinholes)),
        'pinhole_coverage_fraction': float(pinhole_mask.mean()),
        'pinholes': pinholes,
        'al2o3_free_volume_defects': int(len(free_volume_defects)),
        'al2o3_free_volume_defects_detail': free_volume_defects,
        'zno_grain_count': int(n_grains),
        'zno_grain_model': '2d_periodic_voronoi_columns_with_z_wavy_gb_surfaces',
        'zno_gb_fraction_before_blocking': float(gb3d.mean()),
        'gb_block_fraction': float(gb_block_fraction),
        'periodic_axes': ['x', 'y'],
    }
    return grid, nx, ny, nz, al_cells, zno_cells, pinhole_mask, metadata


def summarize(first_times, steps, transmitted, final_z, total_thickness_nm, dx_nm, pinhole_mask):
    ok = transmitted.astype(bool)
    n = int(transmitted.size)
    n_ok = int(ok.sum())
    out = {
        'particles': n,
        'transmitted': n_ok,
        'success_rate': float(n_ok / n) if n else 0.0,
        'pinhole_area_fraction': float(pinhole_mask.mean()),
    }
    if n_ok:
        ft = first_times[ok]
        st = steps[ok]
        out.update({
            'mean_first_passage_s': float(np.mean(ft)),
            'median_first_passage_s': float(np.median(ft)),
            'p10_first_passage_s': float(np.percentile(ft, 10)),
            'p90_first_passage_s': float(np.percentile(ft, 90)),
            'mean_path_steps_success': float(np.mean(st)),
            'median_path_steps_success': float(np.median(st)),
            'D_eff_m2_s': float(((total_thickness_nm * 1e-9) ** 2) / (2.0 * np.mean(ft))),
            'tortuosity_steps_over_thickness': float(np.mean(st) / max(total_thickness_nm / dx_nm, 1.0)),
        })
    else:
        for key in ['mean_first_passage_s', 'median_first_passage_s', 'p10_first_passage_s',
                    'p90_first_passage_s', 'mean_path_steps_success', 'median_path_steps_success',
                    'D_eff_m2_s', 'tortuosity_steps_over_thickness']:
            out[key] = None
    not_ok = ~ok
    if int(not_ok.sum()) > 0:
        depths = final_z[not_ok] * dx_nm
        out['non_transmitted_trapped'] = int(not_ok.sum())
        out['trapped_mean_depth_nm'] = float(np.mean(depths))
        out['trapped_median_depth_nm'] = float(np.median(depths))
    return out


def save_summary_png(out_dir, grid, result, args, pinhole_mask):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    mid_y = grid.shape[2] // 2
    mid_z = grid.shape[0] // 2
    cmap = ListedColormap(['#08306b', '#41b6c4', '#fed976', '#f03b20', '#7a0177', '#333333'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    axes[0, 0].imshow(grid[:, :, mid_y], aspect='auto', origin='upper', cmap=cmap, norm=norm)
    axes[0, 0].set_title('z-x slice at middle y')
    axes[0, 0].set_xlabel('x cell')
    axes[0, 0].set_ylabel('z cell')
    axes[0, 1].imshow(grid[mid_z, :, :].T, origin='lower', cmap=cmap, norm=norm)
    axes[0, 1].set_title('x-y slice at middle z')
    axes[0, 1].set_xlabel('x cell')
    axes[0, 1].set_ylabel('y cell')
    axes[1, 0].imshow(pinhole_mask.T, origin='lower', cmap='gray_r')
    axes[1, 0].set_title('2D surface pinhole mask')
    axes[1, 0].set_xlabel('x cell')
    axes[1, 0].set_ylabel('y cell')
    txt = '\n'.join([
        f"particles = {result['particles']}",
        f"transmitted = {result['transmitted']}",
        f"success_rate = {result['success_rate']:.4f}",
        f"pinhole area frac = {result['pinhole_area_fraction']:.4f}",
        f"mean FPT = {result.get('mean_first_passage_s')} s",
        f"D_eff = {result.get('D_eff_m2_s')} m^2/s",
        f"grid = {grid.shape} (z,x,y)",
        f"elapsed = {result['elapsed_s']:.2f} s",
    ])
    axes[1, 1].text(0.03, 0.97, txt, va='top', family='monospace')
    axes[1, 1].set_axis_off()
    fig.savefig(out_dir / 'bilayer_3d_summary.png', dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='3D on-lattice Al2O3/ZnO bilayer kMC')
    ap.add_argument('--particles', type=int, default=20000)
    ap.add_argument('--max-events', type=int, default=100000)
    ap.add_argument('--max-time-s', type=float, default=1e10)
    ap.add_argument('--width-x-nm', type=float, default=120.0)
    ap.add_argument('--width-y-nm', type=float, default=120.0)
    ap.add_argument('--dx-nm', type=float, default=0.5)
    ap.add_argument('--al-nm', type=float, default=4.5)
    ap.add_argument('--zno-nm', type=float, default=6.0)
    ap.add_argument('--threads-per-block', type=int, default=256)
    ap.add_argument('--out-dir', type=str, default='kmc_bilayer_3d')
    ap.add_argument('--seed', type=int, default=20260529)
    ap.add_argument('--temperature-k', type=float, default=311.15)
    ap.add_argument('--rh', type=float, default=0.90)
    ap.add_argument('--nu', type=float, default=1.0e12)
    ap.add_argument('--e-al-pinhole', type=float, default=0.58)
    ap.add_argument('--e-zno-gb', type=float, default=0.72)
    ap.add_argument('--e-interface', type=float, default=0.62)
    ap.add_argument('--e-surface', type=float, default=0.45)
    ap.add_argument('--gb-block-fraction', type=float, default=0.20)
    ap.add_argument('--skip-plot', action='store_true')
    args = ap.parse_args()

    if not cuda.is_available():
        raise SystemExit('Numba CUDA is not available')

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    grid, nx, ny, nz, al_cells, zno_cells, pinhole_mask, geom_meta = build_bilayer_3d_structure(
        width_x_nm=args.width_x_nm,
        width_y_nm=args.width_y_nm,
        dx_nm=args.dx_nm,
        seed=args.seed,
        al_nm=args.al_nm,
        zno_nm=args.zno_nm,
        gb_block_fraction=args.gb_block_fraction,
    )

    rng = np.random.default_rng(args.seed + 777)
    seeds = rng.integers(1, np.iinfo(np.int64).max, size=args.particles, dtype=np.int64)
    first_times = np.zeros(args.particles, dtype=np.float64)
    steps = np.zeros(args.particles, dtype=np.int32)
    transmitted = np.zeros(args.particles, dtype=np.int32)
    final_z = np.zeros(args.particles, dtype=np.int32)
    final_x = np.zeros(args.particles, dtype=np.int32)
    final_y = np.zeros(args.particles, dtype=np.int32)

    d_grid = cuda.to_device(grid)
    d_mask = cuda.to_device(pinhole_mask)
    d_seeds = cuda.to_device(seeds)
    d_first = cuda.to_device(first_times)
    d_steps = cuda.to_device(steps)
    d_trans = cuda.to_device(transmitted)
    d_fz = cuda.to_device(final_z)
    d_fx = cuda.to_device(final_x)
    d_fy = cuda.to_device(final_y)

    inv_kbt = 1.0 / (KB_EV_K * args.temperature_k)
    blocks = (args.particles + args.threads_per_block - 1) // args.threads_per_block
    bilayer_3d_kernel[blocks, args.threads_per_block](
        d_grid, nx, ny, nz, args.max_events, args.max_time_s, args.nu, inv_kbt,
        args.e_al_pinhole, args.e_zno_gb, args.e_interface, args.e_surface, d_mask,
        d_seeds, d_first, d_steps, d_trans, d_fz, d_fx, d_fy,
    )
    cuda.synchronize()

    first_times = d_first.copy_to_host()
    steps = d_steps.copy_to_host()
    transmitted = d_trans.copy_to_host()
    final_z = d_fz.copy_to_host()
    final_x = d_fx.copy_to_host()
    final_y = d_fy.copy_to_host()
    elapsed = time.time() - t0

    total_thickness_nm = args.al_nm + args.zno_nm
    result = summarize(first_times, steps, transmitted, final_z, total_thickness_nm, args.dx_nm, pinhole_mask)
    result.update({
        'elapsed_s': float(elapsed),
        'model': 'bilayer_3d_on_lattice_v1',
        'environment': {'temperature_c': 38.0, 'temperature_k': args.temperature_k, 'rh': args.rh},
        'total_thickness_nm': float(total_thickness_nm),
        'width_x_nm': float(args.width_x_nm),
        'width_y_nm': float(args.width_y_nm),
        'dx_nm': float(args.dx_nm),
        'max_events': int(args.max_events),
        'barriers_eV': {
            'E_Al_pinhole': args.e_al_pinhole,
            'E_ZnO_GB': args.e_zno_gb,
            'E_interface': args.e_interface,
            'E_surface': args.e_surface,
            'E_Al_matrix': 'infinite (impermeable)',
            'E_ZnO_bulk': 'infinite (impermeable)',
        },
        'geometry_metadata': geom_meta,
    })

    with open(out_dir / 'bilayer_3d_results.csv', 'w', newline='') as f:
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

    with open(out_dir / 'bilayer_3d_results.json', 'w') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    np.savez_compressed(out_dir / 'bilayer_3d_raw.npz', first_times=first_times, steps=steps,
                        transmitted=transmitted, final_z=final_z, final_x=final_x,
                        final_y=final_y, grid=grid, pinhole_mask=pinhole_mask)
    if not args.skip_plot:
        save_summary_png(out_dir, grid, result, args, pinhole_mask)

    print(json.dumps({k: v for k, v in result.items() if k != 'geometry_metadata'}, indent=2, ensure_ascii=False))
    print(f'outputs: {out_dir.resolve()}')


if __name__ == '__main__':
    main()
