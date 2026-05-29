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
    assert period["zno_grain_model"] == "single_truncated_lateral_wurtzite_grain_layer"
    assert period["zno_thickness_to_mean_grain_size_ratio"] < 1.0

    zno = grid[al_cells:height]
    interior = zno[1:-1]
    gb_or_blocked = (interior == model.CELL_ZNO_GB) | (interior == model.CELL_BLOCKED)
    gb_columns = np.where(gb_or_blocked.any(axis=0))[0]
    assert gb_columns.size > 0

    for col in gb_columns:
        ys = np.where(gb_or_blocked[:, col])[0]
        assert ys.size == interior.shape[0]
        assert np.array_equal(ys, np.arange(interior.shape[0]))


if __name__ == "__main__":
    test_sub_10nm_zno_uses_single_truncated_grain_layer()
