from pathlib import Path
import importlib.util

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "kmc_bilayer_all" / "kmc_bilayer_3d.py"


def load_model():
    spec = importlib.util.spec_from_file_location("kmc_bilayer_3d", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_3d_geometry_preserves_10p5nm_thickness_and_periodic_xy_mask():
    model = load_model()
    grid, nx, ny, nz, al_cells, zno_cells, pinhole_mask, metadata = model.build_bilayer_3d_structure(
        width_x_nm=60.0,
        width_y_nm=60.0,
        dx_nm=0.5,
        seed=31,
        al_nm=4.5,
        zno_nm=6.0,
        al_free_volume_count=4,
        gb_block_fraction=0.0,
    )

    assert grid.shape == (21, 120, 120)
    assert nz == 21
    assert al_cells == 9
    assert zno_cells == 12
    assert metadata["total_thickness_nm"] == 10.5
    assert metadata["periodic_axes"] == ["x", "y"]
    assert pinhole_mask.shape == (nx, ny)
    assert pinhole_mask.any()
    assert np.all(grid[0][pinhole_mask.astype(bool)] == model.CELL_AL_PINHOLE)


def test_3d_zno_contains_gb_network_and_interfaces():
    model = load_model()
    grid, nx, ny, nz, al_cells, zno_cells, pinhole_mask, metadata = model.build_bilayer_3d_structure(
        width_x_nm=50.0,
        width_y_nm=50.0,
        dx_nm=0.5,
        seed=43,
        al_nm=4.5,
        zno_nm=6.0,
        al_free_volume_count=2,
        gb_block_fraction=0.10,
    )

    zno = grid[al_cells:nz]
    assert np.all(grid[al_cells - 1][pinhole_mask == 0] == model.CELL_INTERFACE)
    assert np.all(zno[0] == model.CELL_INTERFACE)
    assert np.all(zno[-1] == model.CELL_INTERFACE)
    interior = zno[1:-1]
    assert ((interior == model.CELL_ZNO_GB) | (interior == model.CELL_BLOCKED)).any()
    assert (interior == model.CELL_ZNO_BULK).any()
    assert metadata["zno_grain_model"] == "2d_periodic_voronoi_columns_with_z_wavy_gb_surfaces"
    assert metadata["zno_gb_fraction_before_blocking"] > 0.0


if __name__ == "__main__":
    test_3d_geometry_preserves_10p5nm_thickness_and_periodic_xy_mask()
    test_3d_zno_contains_gb_network_and_interfaces()
