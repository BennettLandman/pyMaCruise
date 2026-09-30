from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_object_image_variance import toads_object_image_variance


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_object_image_variance"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsObjectImageVarianceReferenceHarness.java"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS_SHA256 = "c7deb91ab169c267e063a5c91f4a70c4cb5ed675d28fae7b2094fdfb8da7e319"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
CENTROIDS = np.asarray([[0.5, 0.5, 0.5], [0.31, 0.57, 0.76]], dtype=np.float32)
CLASS_TYPES = ("OBJECT", "OBJECT", "OBJECT")


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).view(np.uint32).astype(np.uint64).tolist()


def test_ordinary_object_image_variance_matches_java_excerpt_float_bits() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert REFERENCE["source_sha256"] == SOURCE_SHA256
    images = _load("images")
    memberships = _load("memberships")
    field = _load("intensity_field")
    assert images.shape == field.shape == (4, 3, 3, 2)
    assert memberships.shape == (4, 3, 3, 3)

    # Confirm each class has real membership on an image face, so a cropped
    # or interior-only scan changes the result instead of accidentally passing.
    boundary = np.zeros(images.shape[:3], dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    assert np.all(np.any(memberships[boundary] > 0.0, axis=0))
    no_boundary_memberships = memberships.copy()
    no_boundary_memberships[boundary] = 0.0

    without_field, without_field_stddev = toads_object_image_variance(
        images, memberships, CENTROIDS, class_types=CLASS_TYPES
    )
    fielded, fielded_stddev = toads_object_image_variance(
        images,
        memberships,
        CENTROIDS,
        class_types=CLASS_TYPES,
        intensity_field=field,
    )
    interior_only, interior_stddev = toads_object_image_variance(
        images, no_boundary_memberships, CENTROIDS, class_types=CLASS_TYPES
    )
    assert _bits(without_field) == REFERENCE["without_field"]["normalized_var_img_bits"]
    assert _bits(without_field_stddev) == REFERENCE["without_field"]["stddev_bits"]
    assert _bits(fielded) == REFERENCE["with_field"]["normalized_var_img_bits"]
    assert _bits(fielded_stddev) == REFERENCE["with_field"]["stddev_bits"]
    assert _bits(np.asarray([np.float32(1.0e10)])) == REFERENCE["without_field"]["inverse_variance_bits"][:1]
    assert without_field_stddev[0] == np.float32(1.0e5)
    assert not np.array_equal(without_field_stddev, interior_stddev)
    assert not np.array_equal(without_field, interior_only)
    assert not np.array_equal(without_field_stddev, fielded_stddev)


@pytest.mark.parametrize("class_types", [("OBJECT", "LESION", "OBJECT"), ("OBJECT", "OUTLIER", "OBJECT"), ("OBJECT",)])
def test_image_variance_adapter_rejects_nonordinary_or_unmatched_class_routes(class_types) -> None:
    images = _load("images")
    memberships = _load("memberships")
    with pytest.raises(ValueError):
        toads_object_image_variance(
            images, memberships, CENTROIDS, class_types=class_types
        )
