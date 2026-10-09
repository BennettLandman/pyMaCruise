from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.toads_crop import crop_toads_images
from pymacruise.toads_post_registration_template_refresh import (
    POST_REGISTRATION_REFRESH_CALL_ORDER,
    toads_post_registration_template_refresh,
)

ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_post_registration_template_refresh"
REFERENCE = ROOT / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsPostRegistrationTemplateRefreshReferenceHarness.java"
HARNESS_SHA256 = "2551f7253307a0e6dd0f0e03c65103e27e3b8b44f3fb0b6738ef3692840e22b2"
FIXTURE_SHA256 = {
    "native_t1_channels_ras.nii.gz": "c6b32c5b8a71d851ddf2db5005d6b8dbc0062c97c01a072d45ff3d7733dd8a36",
    "atlas_template_ras.nii.gz": "2b477b6f71542f2784460e4d5a5c7c10f6b1058d59d8c3ca0025acdaa3efcbeb",
}


def _float_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype=">f4").tobytes(order="C")).hexdigest()


def _label_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype=np.uint8).tobytes(order="C")).hexdigest()


def _assert_float_bits(actual: np.ndarray, expected: list[int]) -> None:
    want = np.asarray([int(value) & 0xFFFFFFFF for value in expected], dtype=np.uint32)
    np.testing.assert_array_equal(np.asarray(actual, dtype=np.float32).view(np.uint32), want)


def test_post_registration_refresh_matches_selected_java_sequence() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert reference["source_hashes"] == {
        "wrapper": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
        "preprocess": "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7",
        "atlas": "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2",
        "rotation": "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0",
        "transform": "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac",
        "image_functions": "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892",
        "numerics": "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb",
    }
    for filename, expected_hash in FIXTURE_SHA256.items():
        assert hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == expected_hash

    image_nii = nib.load(str(ROOT / "native_t1_channels_ras.nii.gz"))
    atlas_nii = nib.load(str(ROOT / "atlas_template_ras.nii.gz"))
    assert aff2axcodes(image_nii.affine) == ("R", "A", "S")
    assert aff2axcodes(atlas_nii.affine) == ("R", "A", "S")
    native = np.asarray(image_nii.dataobj, dtype=np.float32).transpose(3, 0, 1, 2).copy()
    atlas = np.asarray(atlas_nii.dataobj, dtype=np.uint8)
    assert list(native.shape[1:]) == reference["native_shape"]
    assert list(atlas.shape) == reference["template_shape"]

    old_bounds = tuple(reference["old_bounds"])
    padding = int(reference["padding"])
    image_res = tuple(reference["image_resolution"])
    template_res = tuple(reference["template_resolution"])
    signs = tuple(reference["orientation_signs"])
    registered = np.asarray([0.035, -0.018, 0.012, 0.16, -0.11, 0.07], dtype=np.float32)
    _assert_float_bits(registered, reference["registered_transform_bits"])

    # Build the saved pre-registration crop; registration then mutates its shared
    # transform while preserving these already-cropped image samples.
    old_crop = crop_toads_images(
        native,
        old_bounds,
        padding,
        np.zeros(6, dtype=np.float32),
        image_resolution=image_res,
        orientation_signs=signs,
    )
    assert _float_sha256(old_crop.images) == reference["cropped_input_sha256"]
    original_native = native.copy()
    original_atlas = atlas.copy()

    result = toads_post_registration_template_refresh(
        cropped_images=old_crop.images,
        native_shape_xyz=tuple(reference["native_shape"]),
        previous_bounds=old_bounds,
        padding=padding,
        registered_crop_transform=registered,
        topology_template=atlas,
        image_resolution_mm=image_res,
        template_resolution_mm=template_res,
        orientation_signs=signs,
        robust_ranges=np.asarray([[0.0, 1.0], [0.0, 1.0]], dtype=np.float32),
        background_ratio=reference["background_ratio"],
        background_label=reference["background_label"],
    )

    assert result.bounds == tuple(reference["bounds"])
    assert list(result.recropped_images.shape[1:]) == reference["recropped_shape"]
    _assert_float_bits(result.registered_transform, reference["registered_transform_bits"])
    _assert_float_bits(result.restored_transform, reference["restored_transform_bits"])
    _assert_float_bits(result.recropped_transform, reference["recropped_transform_bits"])
    assert _float_sha256(result.native_images) == reference["native_image_sha256"]
    assert _float_sha256(result.recropped_images) == reference["recropped_image_sha256"]
    assert _label_sha256(result.transformed_topology) == reference["transformed_topology_sha256"]
    assert result.call_order == POST_REGISTRATION_REFRESH_CALL_ORDER
    assert list(result.call_order) == reference["call_order"]
    np.testing.assert_array_equal(native, original_native)
    np.testing.assert_array_equal(atlas, original_atlas)
    np.testing.assert_array_equal(registered, np.asarray([0.035, -0.018, 0.012, 0.16, -0.11, 0.07], dtype=np.float32))
