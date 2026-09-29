from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
from skimage.measure import marching_cubes

from pymacruise.toads_dura import (
    apply_dura_output_mask,
    dura_initial_level_sets,
    dura_tgdm_evolve,
    dura_signed_distance_field,
    dynamic_pressure_tgdm_point_update,
    dura_correct_white_matter_topology,
)
from pymacruise.topology import topology_correct
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "p18_distancefield_jar_reference"
    / "distancefield-original-jar.json"
)
DYNAMIC_UPDATE_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "dura_dynamic_update"
    / "java_reference.json"
)
DURA_TGDM_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "dura_tgdm_full_evolution"
    / "java_reference.json"
)
DURA_INITIAL_LEVELSETS_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "dura_initial_levelsets"
    / "java_reference.json"
)


def test_dura_initial_level_sets_match_original_jist_jar() -> None:
    reference = json.loads(DURA_INITIAL_LEVELSETS_FIXTURE.read_text())
    shape = tuple(reference["shape_xyz"])
    x, y, z = np.indices(shape, dtype=np.float64)
    wm = (132.0 + 74.0 * np.sin(0.41 * x + 0.17 * y - 0.2 * z)).astype(np.float32)
    gm = (104.0 + 86.0 * np.cos(0.19 * x - 0.33 * y + 0.27 * z)).astype(np.float32)

    actual = dura_initial_level_sets(
        wm,
        gm,
        membership_level=reference["membership_level"],
    )

    for field_name, actual_field in (
        ("inner_level_set", actual.inner_level_set),
        ("outer_level_set", actual.outer_level_set),
        ("inner_signed_distance", actual.inner_signed_distance),
        ("outer_signed_distance", actual.outer_signed_distance),
    ):
        expected = np.asarray(reference[field_name], dtype=np.float32).reshape(shape)
        np.testing.assert_array_equal(actual_field, expected)


def test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras() -> None:
    shape = (9, 10, 8)
    rng = np.random.default_rng(94021)
    wm = np.linspace(-12.0, 278.0, np.prod(shape), dtype=np.float32)
    rng.shuffle(wm)
    wm = wm.reshape(shape)
    paint = np.zeros(shape, dtype=np.float32)
    # An asymmetric narrow painted path includes the image boundary.
    paint[0:4, 2, 3] = 1.0
    paint[3, 2:7, 3] = 1.0
    paint[3, 6, 3:7] = 1.0
    affine = np.diag([-1.0, -1.2, 1.5, 1.0])
    affine[:3, 3] = (8.0, 11.0, -4.0)
    wm_image = nib.Nifti1Image(wm, affine)
    paint_image = nib.Nifti1Image(paint, affine)

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as caught:
        actual = dura_correct_white_matter_topology(wm_image, paint_image)
    assert len(caught) == 2

    expected = topology_correct(
        wm_image,
        paint_image,
        paint_threshold=0.5,
        lowest_level=0.0,
        highest_level=1.0,
        regularize=0.05,
    ).corrected
    np.testing.assert_array_equal(actual.dataobj, expected.dataobj)
    assert actual.get_data_dtype() == np.dtype(np.float32)
    assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
    assert actual.shape == shape


def test_dura_output_mask_matches_wrapper_threshold_order_and_copies() -> None:
    shape = (2, 2, 2)
    original = np.arange(8, dtype=np.float32).reshape(shape) + 10.0
    white_matter = original + 20.0
    gray_matter = original + 40.0
    level_set = np.asarray([0.0, 1.9, 2.0, 2.1, -1.0, 1.0, 2.0, 3.0], dtype=np.float32).reshape(
        shape
    )
    inputs = (original.copy(), white_matter.copy(), gray_matter.copy(), level_set.copy())

    result = apply_dura_output_mask(
        *inputs,
        mask_translation=2.0,
        level_set_translation=0.25,
    )

    expected_mask = np.asarray(level_set < 2.0, dtype=np.uint8)
    np.testing.assert_array_equal(result.corrected_brain_mask, expected_mask)
    np.testing.assert_array_equal(
        result.stripped_volume, np.where(expected_mask, original, 0.0)
    )
    np.testing.assert_array_equal(
        result.corrected_white_matter, np.where(expected_mask, white_matter, 0.0)
    )
    np.testing.assert_array_equal(
        result.corrected_gray_matter, np.where(expected_mask, gray_matter, 0.0)
    )
    np.testing.assert_array_equal(result.brain_mask_level_set, level_set - np.float32(0.25))
    for actual_input, original_input in zip(inputs, (original, white_matter, gray_matter, level_set)):
        np.testing.assert_array_equal(actual_input, original_input)


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        (np.full((2, 2), 0.0), "3D"),
        (np.full((2, 2, 2), np.nan), "finite"),
    ],
)
def test_dura_output_mask_rejects_invalid_fields(
    replacement: np.ndarray, message: str
) -> None:
    field = np.zeros((2, 2, 2), dtype=np.float32)
    with pytest.raises(ValueError, match=message):
        apply_dura_output_mask(field, field, field, replacement)


def test_dura_tgdm_full_profile_matches_vandy_java_harness() -> None:
    reference = json.loads(DURA_TGDM_FIXTURE.read_text())
    assert "original JIST DistanceField.class" in reference["source"]
    assert reference["distance_field_sha256"] == (
        "92294be2db74351446c3e5898fc35b58497f89ccb38772d4d565885f3f0c35b7"
    )
    shape = tuple(reference["dims"])
    x, y, z = np.indices(shape, dtype=np.float64)
    level_set = (
        np.sqrt((x - 5.0) ** 2 + (y - 4.5) ** 2 + (z - 4.0) ** 2) - 3.2
    ).astype(np.float32)
    reference_level_set = (
        np.sqrt((x - 5.25) ** 2 + (y - 4.25) ** 2 + (z - 4.0) ** 2) - 2.6
    ).astype(np.float32)
    gray_matter_volume = (
        102.0 + 31.0 * np.sin(0.38 * x + 0.20 * y) + 27.0 * np.cos(0.31 * y - 0.25 * z)
    ).astype(np.float32)

    actual = dura_tgdm_evolve(
        level_set,
        gray_matter_volume,
        reference_level_set,
        membership_level=reference["membership_level"],
        alpha=reference["alpha"],
        alignment_force=reference["alignment_force"],
        curvature_force=reference["curvature_force"],
        pressure_force=reference["pressure_force"],
        outer_iterations=reference["outer_iterations"],
        record_passes=True,
    )

    assert actual.completed_steps == 150
    assert actual.completed_outer_passes == 5
    assert actual.landmine_restarts == 0
    assert (
        len(actual.pre_reinitialization_level_sets)
        == len(reference["pre_reinitialization_passes"])
        == 5
    )
    assert len(actual.pass_level_sets) == len(reference["passes"]) == 5
    for actual_pre, expected_pre in zip(
        actual.pre_reinitialization_level_sets, reference["pre_reinitialization_passes"]
    ):
        expected = np.asarray(expected_pre, dtype=np.float32).reshape(shape)
        assert np.array_equal(np.signbit(actual_pre), np.signbit(expected))
        assert float(np.max(np.abs(actual_pre - expected))) <= 0.011
    for actual_pass, expected_pass in zip(actual.pass_level_sets, reference["passes"]):
        expected = np.asarray(expected_pass, dtype=np.float32).reshape(shape)
        assert np.array_equal(np.signbit(actual_pass), np.signbit(expected))
        assert float(np.max(np.abs(actual_pass - expected))) <= 0.0041
    expected_final = np.asarray(reference["final"], dtype=np.float32).reshape(shape)
    np.testing.assert_allclose(actual.level_set, expected_final, rtol=0.0, atol=5e-7)


def test_dura_tgdm_folded_cortex_meets_segmentation_and_surface_tolerances() -> None:
    reference = json.loads(
        (Path(__file__).parent / "data" / "synthetic" / "dura_tgdm_folded_cortex" / "java_reference.json").read_text()
    )
    assert reference["phantom"] == "folded_cortex"
    assert "original JIST DistanceField.class" in reference["source"]
    shape = tuple(reference["dims"])
    x, y, z = np.indices(shape, dtype=np.float64)
    center = tuple((size - 1.0) / 2.0 for size in shape)
    dx, dy, dz = x - center[0], y - center[1], z - center[2]
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    theta = np.arctan2(dy, dx)
    polar = np.arccos(np.divide(dz, radius, out=np.zeros_like(radius), where=radius != 0))
    polar[radius == 0] = 0.0
    folds = 1.15 * np.sin(7.0 * theta) * np.sin(4.0 * polar) + 0.45 * np.cos(
        5.0 * theta + 0.7
    ) * np.sin(3.0 * polar)
    outer = (radius - (8.0 + folds)).astype(np.float32)
    inner = (radius - (4.8 + 0.52 * folds)).astype(np.float32)
    gray_matter = (
        118.0 + 42.0 * np.sin(0.24 * x + 0.16 * y) + 35.0 * np.cos(0.19 * y - 0.21 * z)
    ).astype(np.float32)

    actual = dura_tgdm_evolve(
        outer,
        gray_matter,
        inner,
        membership_level=reference["membership_level"],
        alpha=reference["alpha"],
        alignment_force=reference["alignment_force"],
        curvature_force=reference["curvature_force"],
        pressure_force=reference["pressure_force"],
        outer_iterations=reference["outer_iterations"],
    )
    expected = np.asarray(reference["final"], dtype=np.float32).reshape(shape)

    assert reference["pass_steps"] == [27, 30, 30, 30, 30]
    assert reference["mine_hit_steps"] == [1, 0, 0, 0, 0]
    assert actual.completed_steps == sum(reference["pass_steps"])
    assert actual.landmine_restarts == 1
    assert dice_coefficient(actual.level_set <= 0, expected <= 0) >= 0.97
    actual_vertices, actual_faces, _, _ = marching_cubes(actual.level_set, level=0.0)
    expected_vertices, expected_faces, _, _ = marching_cubes(expected, level=0.0)
    distances = mesh_surface_distance_metrics(
        expected_vertices,
        expected_faces,
        actual_vertices,
        actual_faces,
    )
    assert distances["mean_absolute_distance_mm"] <= 0.01
    assert distances["hausdorff_mm"] <= 1.0


@pytest.mark.parametrize("case_index", range(4))
def test_dura_distance_field_matches_original_jist_jar(case_index: int) -> None:
    reference = json.loads(FIXTURE.read_text())
    case = reference["cases"][case_index]
    actual = dura_signed_distance_field(
        np.asarray(case["input_xyz"], dtype=np.float32),
        max_distance_voxels=case["max_distance_voxels"],
    )

    assert actual.dtype == np.float32
    assert list(actual.shape) == case["shape_xyz"]
    np.testing.assert_array_equal(actual, np.asarray(case["output_xyz"], dtype=np.float32))


@pytest.mark.parametrize(
    ("field", "cap", "message"),
    [
        (np.ones((2, 2), dtype=np.float32), 3.0, "3D"),
        (np.ones((2, 2, 2), dtype=np.float32), 0.0, "positive"),
        (np.full((2, 2, 2), np.nan, dtype=np.float32), 3.0, "finite"),
    ],
)
def test_dura_distance_field_rejects_invalid_inputs(
    field: np.ndarray, cap: float, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        dura_signed_distance_field(field, cap)


@pytest.mark.parametrize("result_index", range(2))
def test_dynamic_pressure_point_update_matches_vandycruise_java(
    result_index: int,
) -> None:
    reference_data = json.loads(DYNAMIC_UPDATE_FIXTURE.read_text())
    i, j, k = np.indices(reference_data["shape"], dtype=np.float64)
    dx, dy, dz = i - 2.0, j - 2.0, k - 2.0
    phi = (0.18 * dx + 0.07 * dy**2 + 0.03 * dy * dz).astype(np.float32)
    reference_level_set = (0.15 * dx + 0.12 * dz + 0.03 * dy**2).astype(np.float32)
    case = reference_data["results"][result_index]
    offset = 0.35 if result_index == 0 else -0.35
    region_forces = (offset + 0.12 * dx - 0.08 * dy + 0.04 * dz).astype(np.float32)
    parameters = reference_data["parameters"]

    actual = dynamic_pressure_tgdm_point_update(
        phi,
        region_forces,
        reference_level_set,
        tuple(case["xyz"]),
        time_step=parameters["time_step"],
        alpha=parameters["alpha"],
        alignment_force=parameters["alignment_force"],
        curvature_force=parameters["curvature_profile_force"],
        pressure_profile_force=parameters["pressure_profile_force"],
    )

    assert actual == pytest.approx(case["updated_phi"], abs=5e-10, rel=0.0)


def test_dynamic_pressure_point_update_rejects_invalid_fields_and_point() -> None:
    field = np.zeros((3, 3, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="identical shapes"):
        dynamic_pressure_tgdm_point_update(
            field,
            field[:2],
            field,
            (1, 1, 1),
            time_step=0.1,
            alpha=0.5,
            alignment_force=1.0,
            curvature_force=-0.5,
            pressure_profile_force=0.5,
        )
    with pytest.raises(ValueError, match="in-bounds"):
        dynamic_pressure_tgdm_point_update(
            field,
            field,
            field,
            (3, 1, 1),
            time_step=0.1,
            alpha=0.5,
            alignment_force=1.0,
            curvature_force=-0.5,
            pressure_profile_force=0.5,
        )
