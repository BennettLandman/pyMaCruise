from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_optimized_centroid_update


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_optimized_centroid" / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"


def _load(directory: str, name: str) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    images = _load("toads_optimized_likelihood", "images")
    field = _load("toads_optimized_likelihood", "field")
    memberships = _load("toads_direct_membership", "initial_memberships")[..., :3]
    return images, field, memberships


def test_optimized_damped_centroids_match_java_power_table_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["power_table_source_sha256"] == POWER_TABLE_SHA256
    assert reference["shape"] == [7, 6, 5, 3]
    centroids = np.asarray(reference["centroids"], dtype=np.float32)
    exponents = np.asarray(reference["exponents"], dtype=np.float32)
    images, field, memberships = _inputs()
    assert images.shape == field.shape == (7, 6, 5, 3)
    assert memberships.shape == (7, 6, 5, 3)
    assert images[3, 3, 2, 0] == centroids[0, 0]  # dist==0 uses source factor 1

    without_field = toads_optimized_centroid_update(
        images, memberships, centroids, exponents, np.float32(0.1)
    )
    with_field = toads_optimized_centroid_update(
        images,
        memberships,
        centroids,
        exponents,
        np.float32(0.1),
        intensity_field=field,
    )
    java_without = np.asarray(reference["without_field"], dtype=np.float32).reshape(centroids.shape)
    java_with = np.asarray(reference["with_field"], dtype=np.float32).reshape(centroids.shape)
    np.testing.assert_array_equal(without_field, java_without)
    np.testing.assert_array_equal(with_field, java_with)
    assert not np.array_equal(without_field, with_field)


def test_optimized_centroid_keeps_previous_value_when_support_is_empty() -> None:
    images = np.full((5, 5, 5, 1), 0.5, dtype=np.float32)
    memberships = np.ones((5, 5, 5, 2), dtype=np.float32)
    memberships[..., 1] = 0.0
    centroids = np.asarray([[0.5, 0.37]], dtype=np.float32)
    exponents = np.ones_like(centroids)
    updated = toads_optimized_centroid_update(
        images, memberships, centroids, exponents, np.float32(0.1)
    )
    assert updated[0, 0] == np.float32(0.5)  # includes exact-zero-distance samples
    assert updated[0, 1] == centroids[0, 1]


@pytest.mark.parametrize(
    ("images", "memberships", "centroids", "exponents", "field"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((1, 1)), np.ones((1, 1)), None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((2, 3, 3, 1)), np.zeros((1, 1)), np.ones((1, 1)), None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((2, 1)), np.ones((1, 1)), None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((1, 1)), np.ones((1, 1)), np.zeros((3, 3, 3))),
    ],
)
def test_optimized_centroid_rejects_invalid_inputs(images, memberships, centroids, exponents, field) -> None:
    with pytest.raises(ValueError):
        toads_optimized_centroid_update(
            images, memberships, centroids, exponents, 0.1, intensity_field=field
        )
