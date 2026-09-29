from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.plus_pipeline import (
    PlusSurfaceResult,
    PlusUncroppedTGDMResult,
    extract_plus_surfaces,
    run_plus_thickness_stage,
)
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness


FIXTURE = Path(__file__).parent / "data/synthetic/plus_surfaces"


def _level_sets() -> PlusUncroppedTGDMResult:
    inner = nib.load(FIXTURE / "inner.nii.gz")
    central = nib.load(FIXTURE / "central.nii.gz")
    outer = nib.load(FIXTURE / "outer.nii.gz")
    return PlusUncroppedTGDMResult(
        initial_level_set=outer,
        tissue_sign=nib.Nifti1Image(
            np.zeros(inner.shape, dtype=np.uint8), inner.affine
        ),
        initial_inner=inner,
        inner=inner,
        central=central,
        outer=outer,
    )


def test_plus_thickness_adapter_routes_volume_and_surface_measurements() -> None:
    level_sets = _level_sets()
    surfaces = extract_plus_surfaces(level_sets)
    actual = run_plus_thickness_stage(level_sets, surfaces)
    direct_volume = levelset_subtraction_thickness(
        level_sets.inner, level_sets.outer, iso_level=0.0
    )
    direct_surface = closest_point_thickness(surfaces.inner, surfaces.outer)

    np.testing.assert_array_equal(actual.volume.image.dataobj, direct_volume.image.dataobj)
    assert nib.aff2axcodes(actual.volume.image.affine) == ("R", "A", "S")
    assert actual.volume.crossing_voxel_count > 0
    assert actual.volume.mean_mm == 2.0
    assert actual.volume.standard_deviation_mm == 0.0
    assert actual.volume.mean_mm == direct_volume.mean_mm

    assert actual.surface.mean_mm == direct_surface.mean_mm
    assert actual.surface.standard_deviation_mm == direct_surface.standard_deviation_mm
    assert actual.surface.nonzero_vertex_count == direct_surface.nonzero_vertex_count
    assert actual.surface.nonzero_vertex_count == len(surfaces.inner.vertices)
    np.testing.assert_array_equal(
        actual.surface.inner_surface.point_data["mean_thickness_mm"],
        direct_surface.inner_surface.point_data["mean_thickness_mm"],
    )
    assert np.isfinite(actual.surface.inner_surface.point_data["mean_thickness_mm"]).all()
    assert actual.surface.inner_surface.coordinate_frame == "RAS"
    assert actual.surface.inner_surface.units == "mm"
