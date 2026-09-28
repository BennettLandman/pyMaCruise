from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_dura import (
    dura_signed_distance_field,
    dynamic_pressure_tgdm_point_update,
)


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
