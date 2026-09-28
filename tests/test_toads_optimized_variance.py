from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_math import toads_optimized_image_variance


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_optimized_variance" / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _load(directory: str, name: str) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def test_optimized_image_variance_matches_java_formula_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["shape"] == [7, 6, 5, 3]
    images = _load("toads_optimized_likelihood", "images")
    field = _load("toads_optimized_likelihood", "field")
    memberships = _load("toads_direct_membership", "initial_memberships")[..., :3]
    centroids = np.asarray(
        [[0.64, 0.62, 0.81], [0.35, 0.55, 0.72], [0.18, 0.70, 0.88]],
        dtype=np.float32,
    )

    without_field = toads_optimized_image_variance(images, memberships, centroids)
    with_field = toads_optimized_image_variance(
        images, memberships, centroids, intensity_field=field
    )
    for actual, key in ((without_field, "without_field"), (with_field, "with_field")):
        np.testing.assert_array_equal(
            actual[0], np.asarray(reference[key]["normalized_inverse_variance"], dtype=np.float32)
        )
        np.testing.assert_array_equal(
            actual[1], np.asarray(reference[key]["image_stddev"], dtype=np.float32)
        )
    assert not np.array_equal(without_field[0], with_field[0])


def test_one_modality_variance_weight_normalizes_to_one() -> None:
    images = np.asarray([[[[0.1], [0.4]], [[0.7], [1.0]]]], dtype=np.float32)
    memberships = np.full((*images.shape[:3], 2), 0.5, dtype=np.float32)
    centroids = np.asarray([[0.2, 0.8]], dtype=np.float32)
    weights, stddev = toads_optimized_image_variance(images, memberships, centroids)
    assert weights.shape == stddev.shape == (1,)
    assert weights[0] == np.float32(1.0)
    assert np.isfinite(stddev).all()


def test_zero_residual_preserves_source_infinite_variance_normalization() -> None:
    images = np.ones((2, 2, 2, 2), dtype=np.float32)
    memberships = np.ones((2, 2, 2, 1), dtype=np.float32)
    centroids = np.ones((2, 1), dtype=np.float32)
    weights, stddev = toads_optimized_image_variance(images, memberships, centroids)
    np.testing.assert_array_equal(weights, np.asarray([0.5, 0.5], dtype=np.float32))
    np.testing.assert_array_equal(stddev, np.full(2, np.float32(1e5)))
