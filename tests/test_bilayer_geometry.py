from pathlib import Path
import importlib.util

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "kmc_bilayer_all" / "kmc_bilayer_corrected.py"


def load_model():
    spec = importlib.util.spec_from_file_location("kmc_bilayer_corrected", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_sub_10nm_zno_uses_single_truncated_grain_layer():
    model = load_model()

    grid, width, height, al_cells, zno_cells, pinhole_mask, metadata = model.build_bilayer_structure(
        width_nm=80.0,
        dx_nm=0.5,
        seed=17,
        periods=1,
        al_nm=4.5,
        zno_nm=6.0,
        grain_min_nm=9.0,
        grain_max_nm=11.0,
        gb_block_fraction=0.0,
    )

    assert zno_cells == 12
    period = metadata["periods"][0]
    assert period["zno_grain_model"] == "single_truncated_lateral_wurtzite_grain_layer_wavy_gb"
    assert period["zno_thickness_to_mean_grain_size_ratio"] < 1.0

    zno = grid[al_cells:height]
    interior = zno[1:-1]
    gb_or_blocked = (interior == model.CELL_ZNO_GB) | (interior == model.CELL_BLOCKED)
    assert gb_or_blocked.any()

    for row in gb_or_blocked:
        assert row.any()

    centers = []
    for yy, row in enumerate(gb_or_blocked):
        xs = np.where(row)[0]
        centers.append(float(xs.mean()))
        if yy > 0:
            prev = np.where(gb_or_blocked[yy - 1])[0]
            assert np.min(np.abs(xs[:, None] - prev[None, :])) <= 3

    assert len(set(round(c, 1) for c in centers)) > 1


def test_al2o3_has_non_through_free_volume_defects():
    model = load_model()

    grid, width, height, al_cells, zno_cells, pinhole_mask, metadata = model.build_bilayer_structure(
        width_nm=120.0,
        dx_nm=0.5,
        seed=29,
        periods=1,
        al_nm=4.5,
        zno_nm=6.0,
        al_free_volume_count=12,
        gb_block_fraction=0.0,
    )

    period = metadata["periods"][0]
    assert period["al2o3_free_volume_defects"] > 0

    al = grid[:al_cells]
    pinhole_cols = np.where(pinhole_mask == 1)[0]
    al_pinhole_cols = np.where((al == model.CELL_AL_PINHOLE).any(axis=0))[0]
    hidden_defect_cols = [col for col in al_pinhole_cols if col not in set(pinhole_cols)]
    assert hidden_defect_cols

    for col in hidden_defect_cols:
        ys = np.where(al[:, col] == model.CELL_AL_PINHOLE)[0]
        assert 0 not in set(ys)
        assert (al_cells - 1) not in set(ys)


if __name__ == "__main__":
    test_sub_10nm_zno_uses_single_truncated_grain_layer()
    test_al2o3_has_non_through_free_volume_defects()
