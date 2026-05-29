#!/usr/bin/env python3
"""Rebuild the readable analysis figure for 3D Al2O3/ZnO bilayer kMC.

Inputs, by default, are expected in an output directory created by
``kmc_bilayer_3d.py``:

- bilayer_3d_raw.npz
- bilayer_3d_results.json

Outputs:

- bilayer_3d_readable_analysis.png
- bilayer_3d_readable_analysis.pdf

The script is intentionally standalone so the published report figure can be
reproduced from the committed raw NPZ + JSON summary without depending on the
one-off notebook/code used during analysis.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, FancyArrowPatch


CELL_NAMES = {
    0: "Al2O3 matrix",
    1: "ZnO bulk/grain",
    2: "blocked",
    3: "Al2O3 pinhole",
    4: "ZnO GB/channel",
    5: "interface",
}

CELL_COLORS = {
    0: "#224c9a",  # Al matrix
    1: "#f2c94c",  # ZnO grain/bulk
    2: "#262626",  # blocked
    3: "#31c8d4",  # pinhole
    4: "#d63b2f",  # GB/channel
    5: "#8e44ad",  # interface
}


def _load_font() -> None:
    """Use an installed CJK font when available; otherwise stay portable."""
    # Keep English labels by default, but this also prevents minus-sign issues.
    mpl.rcParams["axes.unicode_minus"] = False
    candidates = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
    ]
    for fp in candidates:
        if Path(fp).exists():
            try:
                import matplotlib.font_manager as fm

                fm.fontManager.addfont(fp)
                name = fm.FontProperties(fname=fp).get_name()
                mpl.rcParams["font.family"] = name
            except Exception:
                pass
            break


def _first_existing(summary: dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        if key in summary:
            return summary[key]
    return default


def _safe_percentile(x: np.ndarray, q: float) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.percentile(x, q))


def _plot_cell_legend(ax) -> None:
    legend_items = [
        ("Al2O3 dense matrix", CELL_COLORS[0]),
        ("Al2O3 pinhole", CELL_COLORS[3]),
        ("Al2O3/ZnO interface", CELL_COLORS[5]),
        ("ZnO grain/bulk", CELL_COLORS[1]),
        ("ZnO GB/channel", CELL_COLORS[4]),
        ("blocked/dead end", CELL_COLORS[2]),
    ]
    ax.axis("off")
    ax.text(0.0, 1.0, "Cell types", fontsize=12, weight="bold", va="top")
    y = 0.82
    for label, color in legend_items:
        ax.add_patch(Rectangle((0.02, y - 0.035), 0.08, 0.055, color=color,
                               transform=ax.transAxes, clip_on=False))
        ax.text(0.14, y, label, transform=ax.transAxes, va="center", fontsize=10)
        y -= 0.13


def build_figure(raw_npz: Path, summary_json: Path, out_png: Path, out_pdf: Path, dpi: int = 220) -> None:
    _load_font()

    raw = np.load(raw_npz)
    with open(summary_json, "r", encoding="utf-8") as f:
        summary = json.load(f)

    first_times = raw["first_times"]
    steps = raw["steps"]
    transmitted = raw["transmitted"].astype(bool)
    grid = raw["grid"]
    pinhole_mask = raw["pinhole_mask"]

    dx_nm = float(_first_existing(summary, "dx_nm", default=0.5))
    total_thickness_nm = float(_first_existing(summary, "total_thickness_nm", default=grid.shape[0] * dx_nm))
    particles = int(_first_existing(summary, "particles", default=transmitted.size))
    n_trans = int(transmitted.sum())
    ptrans = n_trans / particles if particles else float("nan")
    blocked_pct = 100.0 * (1.0 - ptrans) if math.isfinite(ptrans) else float("nan")

    ok = transmitted
    fpt = first_times[ok]
    path_um = steps[ok].astype(float) * dx_nm / 1000.0
    tortuosity = float(np.mean(path_um) * 1000.0 / total_thickness_nm) if path_um.size else float("nan")
    corr = float(np.corrcoef(path_um, fpt)[0, 1]) if path_um.size > 1 else float("nan")

    # Representative x-y cross-section through the middle lateral plane.
    # grid convention is [z, x, y]; plotting z vertically and x horizontally.
    y_mid = grid.shape[2] // 2
    cross = grid[:, :, y_mid]
    cmap = ListedColormap([CELL_COLORS[i] for i in range(6)])
    norm = BoundaryNorm(np.arange(-0.5, 6.5, 1), cmap.N)

    fig = plt.figure(figsize=(15, 10), facecolor="white")
    gs = fig.add_gridspec(3, 3, width_ratios=[1.3, 1.1, 1.0], height_ratios=[1.0, 1.0, 0.75],
                          wspace=0.36, hspace=0.42)

    fig.suptitle("3D Al2O3/ZnO bilayer kMC: readable transport analysis",
                 fontsize=18, weight="bold", y=0.98)
    fig.text(0.5, 0.945,
             "Rare Al2O3 pinholes + interface lateral search + ZnO boundary/channel tortuosity make successful transmission uncommon.",
             ha="center", fontsize=11, color="#333333")

    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(cross, cmap=cmap, norm=norm, aspect="auto", origin="upper")
    ax0.set_title(f"Middle 3D slice: grid[z, x, y_mid={y_mid}]", fontsize=12, weight="bold")
    ax0.set_xlabel("x grid index")
    ax0.set_ylabel("z thickness index")
    ax0.text(0.02, 0.94, f"thickness = {total_thickness_nm:g} nm\nshape = {grid.shape}",
             transform=ax0.transAxes, va="top", fontsize=9,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.92))

    ax_leg = fig.add_subplot(gs[0, 1])
    _plot_cell_legend(ax_leg)

    ax1 = fig.add_subplot(gs[0, 2])
    ax1.set_title("Transport outcome", fontsize=12, weight="bold")
    vals = [particles - n_trans, n_trans]
    colors = ["#bdbdbd", "#ff8c2a"]
    labels = [f"blocked\n{blocked_pct:.2f}%", f"transmitted\n{ptrans*100:.3f}%"]
    ax1.pie(vals, labels=labels, colors=colors, startangle=90,
            textprops={"fontsize": 10}, wedgeprops={"edgecolor": "white", "linewidth": 1.5})
    ax1.text(0, -1.25, f"{n_trans:,} / {particles:,} transmitted\n≈ 1 in {round(1/ptrans):,} water molecules",
             ha="center", va="top", fontsize=10)

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title("First-passage time distribution", fontsize=12, weight="bold")
    ax2.hist(fpt, bins=45, color="#3f7fc2", alpha=0.85, edgecolor="white")
    ax2.axvline(np.mean(fpt), color="#d62728", lw=2, label=f"mean {np.mean(fpt):.1f} s")
    ax2.axvline(np.median(fpt), color="#111111", lw=2, ls="--", label=f"median {np.median(fpt):.1f} s")
    ax2.set_xlabel("FPT (s), successful trajectories only")
    ax2.set_ylabel("count")
    ax2.legend(frameon=False, fontsize=9)
    ax2.grid(alpha=0.2)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title("Path length vs FPT", fontsize=12, weight="bold")
    # Use deterministic downsampling only if needed; current formal data has 1724 successes.
    rng = np.random.default_rng(20260529)
    if path_um.size > 8000:
        idx = rng.choice(path_um.size, size=8000, replace=False)
    else:
        idx = np.arange(path_um.size)
    ax3.scatter(path_um[idx], fpt[idx], s=10, alpha=0.42, color="#f28e2b", edgecolors="none")
    ax3.set_xlabel("path length (μm) = steps × dx / 1000")
    ax3.set_ylabel("FPT (s)")
    ax3.grid(alpha=0.22)
    ax3.text(0.04, 0.96, f"n = {path_um.size:,}\ncorr = {corr:.3f}",
             transform=ax3.transAxes, va="top", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#dddddd", alpha=0.9))

    ax4 = fig.add_subplot(gs[1, 2])
    ax4.set_title("Path elongation", fontsize=12, weight="bold")
    straight_um = total_thickness_nm / 1000.0
    mean_path_um = float(np.mean(path_um)) if path_um.size else float("nan")
    ax4.barh([1, 0], [straight_um, mean_path_um], color=["#7f7f7f", "#ff8c2a"], height=0.42)
    ax4.set_yticks([1, 0])
    ax4.set_yticklabels([f"straight\n{total_thickness_nm:g} nm", f"tortuous mean\n{mean_path_um:.1f} μm"])
    ax4.set_xlabel("length (μm, log scale)")
    ax4.set_xscale("log")
    ax4.grid(axis="x", alpha=0.25)
    ax4.text(0.98, 0.15, f"tortuosity ≈ {tortuosity:.0f}×",
             transform=ax4.transAxes, ha="right", fontsize=13, weight="bold", color="#c45a00")

    ax5 = fig.add_subplot(gs[2, :])
    ax5.axis("off")
    summary_lines = [
        ("P_trans", f"{ptrans*100:.3f}% ({n_trans:,}/{particles:,})"),
        ("mean FPT", f"{np.mean(fpt):.1f} s"),
        ("FPT p10–p90", f"{_safe_percentile(fpt, 10):.1f}–{_safe_percentile(fpt, 90):.1f} s"),
        ("mean path", f"{mean_path_um:.1f} μm"),
        ("path p10–p90", f"{_safe_percentile(path_um, 10):.1f}–{_safe_percentile(path_um, 90):.1f} μm"),
        ("pinhole area", f"{100*float(_first_existing(summary, 'pinhole_area_fraction', default=float(np.mean(pinhole_mask)))):.4f}%"),
    ]
    xs = [0.02, 0.20, 0.36, 0.54, 0.70, 0.86]
    for x, (k, v) in zip(xs, summary_lines):
        ax5.text(x, 0.78, k, transform=ax5.transAxes, fontsize=9.5, color="#555555")
        ax5.text(x, 0.48, v, transform=ax5.transAxes, fontsize=11.5, weight="bold", color="#111111")

    mechanism = (
        "Mechanism: Most water molecules do not enter because Al2O3 pinholes are rare. "
        "A successful molecule must find a pinhole, search laterally at the Al2O3/ZnO interface, "
        "then follow a long boundary/channel path while avoiding blocked dead ends. "
        "Therefore geometric thickness is only 10.5 nm, but successful paths average ~37.2 μm."
    )
    ax5.text(0.02, 0.10, mechanism, transform=ax5.transAxes, fontsize=10.5,
             bbox=dict(boxstyle="round,pad=0.55", fc="#f7f7f7", ec="#dddddd"), wrap=True)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_dir = Path(__file__).resolve().parent / "kmc_bilayer_3d_formal_200k"
    parser.add_argument("--input-dir", type=Path, default=default_dir,
                        help="Directory containing bilayer_3d_raw.npz and bilayer_3d_results.json")
    parser.add_argument("--raw", type=Path, default=None, help="Override raw NPZ path")
    parser.add_argument("--summary", type=Path, default=None, help="Override JSON summary path")
    parser.add_argument("--out-png", type=Path, default=None, help="Output PNG path")
    parser.add_argument("--out-pdf", type=Path, default=None, help="Output PDF path")
    parser.add_argument("--dpi", type=int, default=220, help="PNG output DPI")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    raw = args.raw or input_dir / "bilayer_3d_raw.npz"
    summary = args.summary or input_dir / "bilayer_3d_results.json"
    out_png = args.out_png or input_dir / "bilayer_3d_readable_analysis.png"
    out_pdf = args.out_pdf or input_dir / "bilayer_3d_readable_analysis.pdf"

    missing = [str(p) for p in [raw, summary] if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required input file(s): " + ", ".join(missing))

    build_figure(raw, summary, out_png, out_pdf, dpi=args.dpi)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()
