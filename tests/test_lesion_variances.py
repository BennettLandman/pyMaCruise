from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_variances


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_variances"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
TYPES = ("object", "lesion", "blackhole", "object")
CLASS_CENTROIDS = np.asarray(
    [[0.15, 0.55, 1.05, 1.65], [0.30, 0.70, 1.15, 1.80]], dtype=np.float32
)
LESION_CENTROID = np.asarray([0.85, 0.95], dtype=np.float32)
BLACK_HOLE_CENTROID = np.asarray([0.45, 0.52], dtype=np.float32)


def _bits(values: np.ndarray) -> list[list[int]] | list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_lesion_class_and_auxiliary_variances_match_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "variance_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 11)
    labels = frames[..., 0].astype(np.int16)
    memberships = frames[..., 1:5]
    values = frames[..., 5:7]
    fields = frames[..., 7:9]
    lesion = frames[..., 9]
    black_hole = frames[..., 10]
    lesion_types = labels == 1
    black_types = labels == 2
    assert np.any(lesion_types & (memberships[..., 1] < lesion) & (lesion > black_hole))
    assert np.any(lesion_types & (memberships[..., 1] < lesion) & (lesion < black_hole))
    assert np.any(black_types & (memberships[..., 2] < black_hole) & (lesion >= black_hole))
    assert np.any(black_types & (memberships[..., 2] < black_hole) & (lesion < black_hole))
    assert np.any(black_types & (lesion == black_hole))
    assert np.any((lesion == 0.0) & (black_hole == 0.0))

    no_field = compute_lesion_toads_variances(
        values, memberships, labels, CLASS_CENTROIDS, lesion, black_hole,
        LESION_CENTROID, BLACK_HOLE_CENTROID, TYPES,
    )
    with_field = compute_lesion_toads_variances(
        values, memberships, labels, CLASS_CENTROIDS, lesion, black_hole,
        LESION_CENTROID, BLACK_HOLE_CENTROID, TYPES, fields=fields,
    )
    for actual, expected in ((no_field, REFERENCE["no_field"]), (with_field, REFERENCE["with_field"])):
        assert _bits(actual.class_standard_deviation) == expected["class_bits"]
        assert _bits(actual.lesion_standard_deviation) == expected["lesion_bits"]
        assert _bits(actual.black_hole_standard_deviation) == expected["black_hole_bits"]
    np.testing.assert_array_equal(no_field.class_standard_deviation[:, 3], np.zeros(2, dtype=np.float32))


def test_lesion_variances_reject_malformed_inputs() -> None:
    volume = np.zeros((3, 3, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="object_types"):
        compute_lesion_toads_variances(
            volume[..., None], np.ones((*volume.shape, 2)), np.zeros_like(volume, dtype=np.int16),
            np.zeros((1, 2)), volume, volume, np.zeros(1), np.zeros(1), ("object",),
        )
