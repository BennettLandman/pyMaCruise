from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_preprocessing_transform import transform_lesion_toads_image


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_transform_new_image"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsTransformNewImageReferenceHarness.java"
REFERENCE = FIXTURE_DIR / "java_reference.json"
SOURCE_SHA256 = "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"
IMAGE_FUNCTIONS_SHA256 = "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892"
ROTATION_MATRIX_SHA256 = "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0"
HARNESS_SHA256 = "6676ff3cc7393e040e4affa9bd3e9e009701395e092fdfeccf605f747e6ebca0"
REFERENCE_SHA256 = "ab8eaedcf39dbec3e686509827eeb5773a91a461320f5c903a39f144d49e6b00"


def test_transform_new_image_matches_selected_java_float_bits() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadPreprocess.transformNewImage"
    assert ref["source_sha256"] == SOURCE_SHA256
    assert ref["image_functions_sha256"] == IMAGE_FUNCTIONS_SHA256
    assert ref["rotation_matrix_sha256"] == ROTATION_MATRIX_SHA256

    image = load_nifti_ras(FIXTURE_DIR / "template_image.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.9, 1.0, 1.1), atol=1e-6)
    np.testing.assert_array_equal(image.affine[:3, 3], (12.0, -7.0, 4.0))
    source = np.asarray(image.dataobj, dtype=np.float32)
    assert source.shape == tuple(ref["template_shape"])
    before = source.copy()

    output = transform_lesion_toads_image(
        source,
        (0.02, -0.03, 0.01, 0.3, -0.15, 0.4),
        (7, 8, 9),
        (1, 5, 1, 6, 1, 7),
        image_resolution=(1.1, 1.3, 0.9),
        template_resolution=(0.9, 1.0, 1.1),
        padding=1,
        background=-2.5,
    )
    assert output.shape == tuple(ref["output_shape"])
    expected_bits = np.asarray(ref["output_bits"], dtype=np.uint32).reshape(output.shape)
    np.testing.assert_array_equal(output.view(np.uint32), expected_bits)
    np.testing.assert_array_equal(source, before)
    background = np.float32(-2.5)
    # Crop padding is initialized to background before the source sampling loop.
    assert np.all(output[0] == background) and np.all(output[-1] == background)
    assert np.all(output[:, 0] == background) and np.all(output[:, -1] == background)
    assert np.all(output[:, :, 0] == background) and np.all(output[:, :, -1] == background)
    # The baseline crop's entire six-face pad remains at the background.
    # A separate translated case below exposes the distinct interpolation
    # out-of-domain background path.
    out_of_domain = transform_lesion_toads_image(
        source,
        (0.02, -0.03, 0.01, 2.0, -0.15, 0.4),
        (7, 8, 9),
        (1, 5, 1, 6, 1, 7),
        image_resolution=(1.1, 1.3, 0.9),
        template_resolution=(0.9, 1.0, 1.1),
        padding=1,
        background=-2.5,
    )
    out_of_domain_bits = np.asarray(ref["out_of_domain_bits"], dtype=np.uint32).reshape(
        out_of_domain.shape
    )
    np.testing.assert_array_equal(out_of_domain.view(np.uint32), out_of_domain_bits)
    interior = output[1:-1, 1:-1, 1:-1]
    assert np.count_nonzero(interior != background) > 0
    translated_interior = out_of_domain[1:-1, 1:-1, 1:-1]
    assert np.count_nonzero(translated_interior == background) > 0
    assert np.count_nonzero(translated_interior != background) > 0

    # The stored NIfTI remains RAS; this separately exercises the selected
    # source's signed orientation factor for one sampling axis.
    negative_x = transform_lesion_toads_image(
        source,
        (0.02, -0.03, 0.01, 0.3, -0.15, 0.4),
        (7, 8, 9),
        (1, 5, 1, 6, 1, 7),
        image_resolution=(1.1, 1.3, 0.9),
        template_resolution=(0.9, 1.0, 1.1),
        orientation_signs=(-1.0, 1.0, 1.0),
        padding=1,
        background=-2.5,
    )
    negative_expected = np.asarray(ref["negative_x_orientation_bits"], dtype=np.uint32).reshape(
        negative_x.shape
    )
    np.testing.assert_array_equal(negative_x.view(np.uint32), negative_expected)
    assert not np.array_equal(output.view(np.uint32), negative_x.view(np.uint32))


def test_transform_new_image_rejects_bad_crop_bounds() -> None:
    with np.testing.assert_raises_regex(ValueError, "inside image_shape"):
        transform_lesion_toads_image(
            np.ones((4, 4, 4), dtype=np.float32),
            np.zeros(6, dtype=np.float32),
            (3, 3, 3),
            (0, 3, 0, 2, 0, 2),
            image_resolution=(1, 1, 1),
            template_resolution=(1, 1, 1),
        )
