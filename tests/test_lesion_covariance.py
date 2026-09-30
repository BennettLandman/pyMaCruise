"""Selected-Java comparison for the dormant LesionTOADS covariance method."""

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_covariance import compute_lesion_toads_covariances


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_covariance"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
TYPES = ("object", "object", "lesion", "blackhole", "outlier", "object")
LESION_CENTROID = np.asarray([0.25, 0.5625, 0.8125], dtype=np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _expected(values: list[list[int]]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32)


def test_covariances_match_selected_java_method_bitwise() -> None:
    image = nib.load(str(ROOT / "covariance_inputs.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    packed = np.asarray(image.dataobj, dtype=np.float32)
    assert packed.shape == (5, 4, 3, 14)
    labels = packed[..., 0].astype(np.int16)
    memberships = packed[..., 1:7].copy()
    images = packed[..., 7:10].copy()
    fields = packed[..., 10:13].copy()
    lesion_membership = packed[..., 13].copy()
    centroids = np.asarray(
        [[((c * 3 + k * 2 + 1) % 13) * 0.0625 for k in range(6)] for c in range(3)],
        dtype=np.float32,
    )
    reference = json.loads((ROOT / "java_reference.json").read_text())
    fixture_path = ROOT / "covariance_inputs.nii.gz"
    digest = hashlib.sha256(fixture_path.read_bytes()).hexdigest()
    assert reference["source_sha256"] == SOURCE_SHA
    harness = Path(__file__).parents[1] / "tests" / "java_reference" / "LesionCovarianceReferenceHarness.java"
    assert reference["harness_sha256"] == hashlib.sha256(harness.read_bytes()).hexdigest()
    assert reference["input_sha256"] == digest

    originals = tuple(array.copy() for array in (images, memberships, labels, fields, lesion_membership))
    for mode, field_values in (("no_field", None), ("with_field", fields)):
        result = compute_lesion_toads_covariances(
            images,
            memberships,
            labels,
            centroids,
            lesion_membership,
            LESION_CENTROID,
            TYPES,
            fields=field_values,
        )
        expected = reference[mode]
        np.testing.assert_array_equal(
            result.class_covariance.view(np.int32),
            _expected(expected["class_bits"]).reshape(6, 3, 3).view(np.int32),
        )
        np.testing.assert_array_equal(
            result.lesion_covariance.view(np.int32),
            _expected(expected["lesion_bits"]).reshape(3, 3).view(np.int32),
        )
        np.testing.assert_array_equal(result.class_covariance, result.class_covariance.transpose(0, 2, 1))
        np.testing.assert_array_equal(result.lesion_covariance, result.lesion_covariance.T)
    assert not np.array_equal(
        _expected(reference["no_field"]["lesion_bits"]),
        _expected(reference["with_field"]["lesion_bits"]),
    )
    assert not np.any(_expected(reference["no_field"]["class_bits"]).reshape(6, 3, 3)[5])
    for original, current in zip(originals, (images, memberships, labels, fields, lesion_membership)):
        np.testing.assert_array_equal(current, original)


def test_covariance_rejects_nonmatching_fixed_input_shapes() -> None:
    with np.testing.assert_raises_regex(ValueError, "fields must match"):
        compute_lesion_toads_covariances(
            np.ones((3, 3, 3, 2), dtype=np.float32),
            np.ones((3, 3, 3, 1), dtype=np.float32),
            np.zeros((3, 3, 3), dtype=np.int8),
            np.zeros((2, 1), dtype=np.float32),
            np.zeros((3, 3, 3), dtype=np.float32),
            np.zeros(2, dtype=np.float32),
            ("OBJECT",),
            fields=np.ones((3, 3, 3, 1), dtype=np.float32),
        )


def test_lesion_ties_and_absent_lesion_winners_leave_auxiliary_matrix_zero() -> None:
    image = nib.load(str(ROOT / "covariance_inputs.nii.gz"))
    packed = np.asarray(image.dataobj, dtype=np.float32)
    labels = packed[..., 0].astype(np.int16)
    memberships = packed[..., 1:7].copy()
    lesion_membership = packed[..., 13].copy()
    lesion_membership[labels == 2] = memberships[..., 2][labels == 2] * np.float32(0.5)
    centroids = np.asarray(
        [[((c * 3 + k * 2 + 1) % 13) * 0.0625 for k in range(6)] for c in range(3)],
        dtype=np.float32,
    )
    result = compute_lesion_toads_covariances(
        packed[..., 7:10],
        memberships,
        labels,
        centroids,
        lesion_membership,
        LESION_CENTROID,
        TYPES,
    )
    np.testing.assert_array_equal(result.lesion_covariance, np.zeros((3, 3), dtype=np.float32))
