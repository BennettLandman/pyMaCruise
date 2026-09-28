from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_optimized_class_standard_deviation


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_optimized_class_variance" / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _load(directory: str, name: str) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def test_optimized_class_standard_deviations_match_java_formula_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["shape"] == [7, 6, 5, 3]
    images = _load("toads_optimized_likelihood", "images")
    field = _load("toads_optimized_likelihood", "field")
    memberships = _load("toads_direct_membership", "initial_memberships")[..., :3]
    labels = _load("toads_optimized_class_variance", "segmentation").astype(np.int16)
    centroids = np.asarray(
        [[0.64, 0.62, 0.81], [0.35, 0.55, 0.72], [0.18, 0.70, 0.88]],
        dtype=np.float32,
    )
    without_field = toads_optimized_class_standard_deviation(
        images, memberships, centroids, labels
    )
    with_field = toads_optimized_class_standard_deviation(
        images, memberships, centroids, labels, intensity_field=field
    )
    np.testing.assert_array_equal(
        without_field,
        np.asarray(reference["without_field"], dtype=np.float32).reshape(3, 3),
    )
    np.testing.assert_array_equal(
        with_field,
        np.asarray(reference["with_field"], dtype=np.float32).reshape(3, 3),
    )
    assert not np.array_equal(without_field, with_field)


def test_empty_assigned_class_has_zero_standard_deviation() -> None:
    images = np.arange(27, dtype=np.float32).reshape(3, 3, 3, 1)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    centroids = np.asarray([[0.0, 1.0]], dtype=np.float32)
    labels = np.zeros((3, 3, 3), dtype=np.int16)
    result = toads_optimized_class_standard_deviation(
        images, memberships, centroids, labels
    )
    assert result[0, 0] > 0.0
    assert result[0, 1] == np.float32(0.0)


@pytest.mark.parametrize(
    ("images", "memberships", "centroids", "labels"),
    [
        (np.zeros((3, 3, 3)), np.zeros((3, 3, 3, 1)), np.zeros((1, 1)), np.zeros((3, 3, 3), dtype=np.int16)),
        (np.zeros((3, 3, 3, 1)), np.zeros((2, 3, 3, 1)), np.zeros((1, 1)), np.zeros((3, 3, 3), dtype=np.int16)),
        (np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.zeros((2, 1)), np.zeros((3, 3, 3), dtype=np.int16)),
        (np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.zeros((1, 1)), np.zeros((3, 3, 3), dtype=np.float32)),
    ],
)
def test_optimized_class_standard_deviation_rejects_invalid_inputs(
    images, memberships, centroids, labels
) -> None:
    with pytest.raises(ValueError):
        toads_optimized_class_standard_deviation(images, memberships, centroids, labels)
