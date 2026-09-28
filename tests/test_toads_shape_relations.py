from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_shape_relation_factors


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_shape_relations"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
NUMERICS_SHA256 = "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
PRIOR_FACTOR = np.float32(0.75)


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def _expected_inputs() -> tuple[np.ndarray, np.ndarray]:
    x, y, z, class_index = np.indices((7, 6, 5, 4), dtype=np.int32)
    shape = (((x + 2 * y + 3 * z + 4 * class_index) % 11) / 10.0).astype(np.float32)
    labels = np.full(shape.shape, -1, dtype=np.int16)
    labels[..., 0] = ((x[..., 0] + y[..., 0] + z[..., 0]) % 4).astype(np.int16)
    direct2 = ((x[..., 0] + y[..., 0] + z[..., 0]) % 3) == 0
    labels[..., 1][direct2] = (labels[..., 0][direct2] + 1) % 4
    labels[..., 2] = ((2 * x[..., 0] + y[..., 0]) % 4).astype(np.int16)
    labels[..., 3] = ((x[..., 0] + 2 * y[..., 0] + z[..., 0] + 1) % 4).astype(np.int16)
    return shape, labels


def test_shape_relation_factors_match_java_source_formula() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["numerics_source_sha256"] == NUMERICS_SHA256
    assert reference["shape"] == [7, 6, 5, 4]
    assert np.float32(reference["prior_factor"]) == PRIOR_FACTOR
    shapes, labels = _load("shape_priors"), _load("closest_labels")
    expected_shapes, expected_labels = _expected_inputs()
    np.testing.assert_array_equal(shapes, expected_shapes)
    np.testing.assert_array_equal(labels, expected_labels)

    factors, suppressed = toads_shape_relation_factors(shapes, labels, PRIOR_FACTOR)
    java_factors = np.asarray(reference["factors"], dtype=np.float32).reshape(factors.shape)
    java_suppressed = np.asarray(reference["suppressed"], dtype=np.bool_).reshape(suppressed.shape)
    np.testing.assert_array_equal(factors, java_factors)
    np.testing.assert_array_equal(suppressed, java_suppressed)
    assert np.any(suppressed)
    assert np.any(factors == np.float32(0.25))
    assert np.any(factors == np.float32(1.0))


def test_relation_precedence_and_zero_shape_rules() -> None:
    shapes = np.zeros((1, 1, 1, 3), dtype=np.float32)
    labels = np.full((1, 1, 1, 4), -1, dtype=np.int16)
    # Direct has precedence over a competing neighbor slot.
    labels[0, 0, 0, 0] = 0
    labels[0, 0, 0, 2] = 0
    labels[0, 0, 0, 3] = 1
    shapes[0, 0, 0, 1] = np.float32(0.1)
    shapes[0, 0, 0, 2] = np.float32(0.9)
    factors, suppressed = toads_shape_relation_factors(shapes, labels, PRIOR_FACTOR)
    assert factors[0, 0, 0, 0] == np.float32(1.0)
    assert not suppressed[0, 0, 0, 0]
    assert factors[0, 0, 0, 1] == np.float32(0.25)
    assert not suppressed[0, 0, 0, 1]
    assert factors[0, 0, 0, 2] == np.float32(1.0)
    assert not suppressed[0, 0, 0, 2]

    labels.fill(-1)
    factors, suppressed = toads_shape_relation_factors(shapes, labels, PRIOR_FACTOR)
    assert factors[0, 0, 0, 0] == 0.0
    assert suppressed[0, 0, 0, 0]
    assert factors[0, 0, 0, 1] == np.float32(0.15)
    assert not suppressed[0, 0, 0, 1]
    assert factors[0, 0, 0, 2] == np.float32(1.0)
    assert not suppressed[0, 0, 0, 2]


@pytest.mark.parametrize(
    ("shapes", "labels", "prior_factor"),
    [
        (np.zeros((2, 2, 2), dtype=np.float32), np.zeros((2, 2, 2, 4), dtype=np.int8), 0.5),
        (np.zeros((2, 2, 2, 2), dtype=np.float32), np.zeros((2, 2, 2, 3), dtype=np.int8), 0.5),
        (np.zeros((2, 2, 2, 2), dtype=np.float32), np.zeros((2, 2, 2, 4), dtype=np.float32), 0.5),
        (np.full((2, 2, 2, 2), np.nan, dtype=np.float32), np.zeros((2, 2, 2, 4), dtype=np.int8), 0.5),
        (np.zeros((2, 2, 2, 2), dtype=np.float32), np.full((2, 2, 2, 4), 2, dtype=np.int8), 0.5),
        (np.zeros((2, 2, 2, 2), dtype=np.float32), np.zeros((2, 2, 2, 4), dtype=np.int8), np.nan),
    ],
)
def test_shape_relation_factors_reject_invalid_inputs(shapes, labels, prior_factor) -> None:
    with pytest.raises(ValueError):
        toads_shape_relation_factors(shapes, labels, prior_factor)
