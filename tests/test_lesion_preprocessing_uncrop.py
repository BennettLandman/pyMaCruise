from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_preprocessing_uncrop import (
    LESION_PREPROCESS_SOURCE_SHA256,
    LESION_UNCROP_METHOD_NORMALIZED_SHA256,
    uncrop_lesion_toads_preprocessing_inputs,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_preprocessing_uncrop"
REFERENCE = Path(__file__).parent / "data" / "synthetic" / "toads_crop_geometry" / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsCropGeometryReferenceHarness.java"
REFERENCE_SHA256 = "27de97cd1cee967dadd621a26d760d2a7cd09ad9aa0cdb9680644367c6dbc29a"
HARNESS_SHA256 = "73434da80a7c4d05aff37738f410bd40a617cea1b588253388dc3133e120c4c3"


def test_lesion_uncrop_source_provenance_matches_shared_kernel() -> None:
    provenance = json.loads((ROOT / "source_provenance.json").read_text())
    assert provenance["source_file_sha256"] == LESION_PREPROCESS_SOURCE_SHA256
    assert provenance["normalized_method_body_sha256"] == LESION_UNCROP_METHOD_NORMALIZED_SHA256
    assert (
        provenance["normalized_method_body_sha256"]
        == provenance["matching_kernel_normalized_method_body_sha256"]
    )
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    assert provenance["unmodeled_lifecycle_fields"] == [
        "xCi = nix/2, yCi = niy/2, zCi = niz/2 internal-center updates",
        "cropped = false internal-state update",
    ]
    assert provenance["source_operation_sequence"] == [
        "return immediately if cropped is false",
        "reset xCi/yCi/zCi to original-grid centers",
        "subtract rotation-scaled crop-center offset from transform translation",
        "allocate full-size output and initialize it to zero",
        "copy the clipped inclusive padded region from cropped images",
        "replace images and clear cropped flag",
    ]


def test_lesion_uncrop_reuses_java_validated_ras_crop_geometry_result() -> None:
    case = json.loads(REFERENCE.read_text())
    cropped = np.asarray(case["cropped_image"], dtype=np.float32)[None, ...]
    transform = np.asarray(case["cropped_transform"], dtype=np.float32)
    cropped_before = cropped.copy()
    transform_before = transform.copy()
    assert case["image_orientation_signs"] == [-1.0, -1.0, 1.0]
    assert case["image_resolution"] != case["template_resolution"]
    assert np.any(transform[:3] != 0) and np.any(transform[3:] != 0)

    restored, restored_transform = uncrop_lesion_toads_preprocessing_inputs(
        cropped,
        tuple(case["dims_image"]),
        tuple(case["bounds_inclusive"]),
        int(case["padding"]),
        transform,
        image_resolution_mm=tuple(case["image_resolution"]),
        orientation_signs=tuple(case["image_orientation_signs"]),
    )

    np.testing.assert_array_equal(
        restored[0], np.asarray(case["uncropped_image"], dtype=np.float32)
    )
    expected_transform = np.asarray(case["uncropped_transform"], dtype=np.float32)
    np.testing.assert_array_equal(restored_transform.view(np.int32), expected_transform.view(np.int32))
    assert np.count_nonzero(restored[:, 0]) == 0
    assert np.count_nonzero(restored[:, -1]) == 0
    np.testing.assert_array_equal(cropped, cropped_before)
    np.testing.assert_array_equal(transform, transform_before)
