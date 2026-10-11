from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.ace import _distance_and_characteristic_fmm, _tgdm_signed_distance_fmm


ROOT = Path(__file__).parent
JIST_DISTANCE_FIXTURE = (
    ROOT / "data" / "synthetic" / "p18_distancefield_jar_reference"
    / "distancefield-original-jar.json"
)


def test_tgdm_fmm_matches_all_original_jist_distancefield_cases_exactly() -> None:
    reference = json.loads(JIST_DISTANCE_FIXTURE.read_text())
    assert len(reference["cases"]) == 4
    for case in reference["cases"]:
        input_field = np.asarray(case["input_xyz"], dtype=np.float32)
        expected = np.asarray(case["output_xyz"], dtype=np.float32)
        actual = _tgdm_signed_distance_fmm(
            input_field, max_distance=case["max_distance_voxels"]
        )
        assert actual.dtype == np.float32
        np.testing.assert_array_equal(actual, expected, err_msg=case["name"])


def test_tgdm_fmm_finite_band_equals_full_field_clipping() -> None:
    x = np.arange(24, dtype=np.float32)[:, None, None]
    y = np.arange(21, dtype=np.float32)[None, :, None]
    z = np.arange(18, dtype=np.float32)[None, None, :]
    field = np.sqrt((x - 11.2) ** 2 + (y - 9.4) ** 2 + (z - 8.1) ** 2) - 5.7

    full = _tgdm_signed_distance_fmm(field, max_distance=30.0)
    limited = _tgdm_signed_distance_fmm(field, max_distance=3.0)

    np.testing.assert_array_equal(limited, np.clip(full, -3.0, 3.0))


def test_tgdm_fmm_matches_python_heap_on_folded_synthetic_field() -> None:
    shape = (48, 47, 43)
    x = np.linspace(-1.0, 1.0, shape[0], dtype=np.float32)[:, None, None]
    y = np.linspace(-1.0, 1.0, shape[1], dtype=np.float32)[None, :, None]
    z = np.linspace(-1.0, 1.0, shape[2], dtype=np.float32)[None, None, :]
    radius = np.sqrt((x / 0.8) ** 2 + (y / 0.75) ** 2 + (z / 0.7) ** 2)
    folds = 0.07 * np.sin(7 * x + 2 * y) * np.sin(6 * z - 2 * y)
    folds += 0.03 * np.cos(9 * y + 3 * z)
    field = np.asarray((radius - 0.72 - folds) * 40.0, dtype=np.float32)

    reference, _ = _distance_and_characteristic_fmm(
        field, np.ones(shape, dtype=np.float32), max_distance=6.2
    )
    actual = _tgdm_signed_distance_fmm(field, max_distance=6.2)

    np.testing.assert_array_equal(actual, reference)
