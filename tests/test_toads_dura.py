from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_dura import dura_signed_distance_field


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "p18_distancefield_jar_reference"
    / "distancefield-original-jar.json"
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
