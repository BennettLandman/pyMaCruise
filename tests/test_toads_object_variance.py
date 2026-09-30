from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_object_variance import toads_object_assigned_class_standard_deviation


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_object_variance"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsObjectVarianceReferenceHarness.java"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS_SHA256 = "4a635e69269be3ceebc24dc7c459bf49b90874291d90889533d5897b23671e44"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
CENTROIDS = np.asarray(
    [[0.31, 0.73, 1.1, 1.37], [0.43, 0.89, 1.2, 1.49]], dtype=np.float32
)


def _load(name: str, *, integer: bool = False) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.int16 if integer else np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).view(np.uint32).astype(np.uint64).tolist()


def test_ordinary_object_assigned_class_variances_match_java_excerpt_bits() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert REFERENCE["source_sha256"] == SOURCE_SHA256
    assert np.asarray(_bits(CENTROIDS)).reshape(-1).tolist() == REFERENCE["centroid_bits"]
    images = _load("images")
    field = _load("intensity_field")
    memberships = _load("memberships")
    segmentation = _load("segmentation", integer=True)
    assert images.shape == field.shape == (4, 3, 3, 2)
    assert memberships.shape == (4, 3, 3, 4)
    assert nib.aff2axcodes(nib.load(str(ROOT / "images.nii.gz")).affine) == ("R", "A", "S")
    assert set(np.unique(segmentation)) == {0, 1, 2}  # class 2 is assigned only at boundary corners
    assert np.all(memberships[..., 1] == 0.0)  # assigned class uses the source 1e-30 weight floor
    assert not np.any(segmentation == 3)  # empty assigned class returns zero

    without_field = toads_object_assigned_class_standard_deviation(
        images, memberships, CENTROIDS, segmentation
    )
    with_field = toads_object_assigned_class_standard_deviation(
        images, memberships, CENTROIDS, segmentation, intensity_field=field
    )
    assert np.asarray(_bits(without_field)).reshape(-1).tolist() == REFERENCE["without_field_bits"]
    assert np.asarray(_bits(with_field)).reshape(-1).tolist() == REFERENCE["with_field_bits"]
    assert np.array_equal(without_field[:, 3], np.zeros(2, dtype=np.float32))
    assert np.all(without_field[:, 2] >= 0.0)  # boundary-only class participates in variance
    assert np.all(without_field[:, 1] > 0.0)  # zero memberships remain weighted by 1e-30f
    assert not np.array_equal(without_field[:, :2], with_field[:, :2])


@pytest.mark.parametrize(
    ("images", "memberships", "centroids", "segmentation"),
    [
        (np.zeros((2, 2, 2)), np.zeros((2, 2, 2, 1)), np.zeros((1, 1)), np.zeros((2, 2, 2), dtype=np.int16)),
        (np.zeros((2, 2, 2, 1)), np.zeros((2, 2, 2, 1)), np.zeros((2, 1)), np.zeros((2, 2, 2), dtype=np.int16)),
        (np.zeros((2, 2, 2, 1)), np.zeros((2, 2, 2, 1)), np.zeros((1, 1)), np.zeros((2, 2, 2), dtype=np.float32)),
    ],
)
def test_ordinary_object_variance_rejects_invalid_shapes_and_labels(
    images, memberships, centroids, segmentation
) -> None:
    with pytest.raises(ValueError):
        toads_object_assigned_class_standard_deviation(
            images, memberships, centroids, segmentation
        )
