from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

from pymacruise.toads_dura import dura_tgdm_evolve
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics


FIXTURE = Path(__file__).parent / "data/synthetic/dura_tgdm_folded_cortex_41"


def _load(path: Path) -> np.ndarray:
    image = nib.load(path)
    assert image.shape == (41, 41, 41)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj, dtype=np.float32)


def test_dura_tgdm_matches_folded_41_java_segmentation_and_surface_gates() -> None:
    outer = _load(FIXTURE / "inputs/outer_level_set.nii.gz")
    inner = _load(FIXTURE / "inputs/inner_level_set.nii.gz")
    gray = _load(FIXTURE / "inputs/gray_matter.nii.gz")
    result = dura_tgdm_evolve(
        outer,
        gray,
        inner,
        membership_level=0.4,
        alpha=0.5,
        alignment_force=1.0,
        curvature_force=0.5,
        pressure_force=0.5,
        outer_iterations=5,
        record_passes=True,
    )

    assert result.completed_steps == 150
    assert result.completed_outer_passes == 5
    assert result.landmine_restarts == 0
    for index, actual in enumerate(result.pass_level_sets, start=1):
        expected = _load(FIXTURE / f"java_reference/pass_{index:02d}.nii.gz")
        assert np.array_equal(np.signbit(actual), np.signbit(expected))

    expected = _load(FIXTURE / "java_reference/final.nii.gz")
    actual = np.asarray(result.level_set, dtype=np.float32)
    assert dice_coefficient(actual <= 0, expected <= 0) >= 0.97
    actual_vertices, actual_faces, _, _ = marching_cubes(actual, level=0.0)
    expected_vertices, expected_faces, _, _ = marching_cubes(expected, level=0.0)
    surface = mesh_surface_distance_metrics(
        actual_vertices, actual_faces, expected_vertices, expected_faces
    )
    assert surface["mean_absolute_distance_mm"] <= 0.01
    assert surface["hausdorff_mm"] <= 1.0
