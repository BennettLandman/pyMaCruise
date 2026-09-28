from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_math import toads_optimized_tissue_data_term


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_optimized_likelihood"
REFERENCE = ROOT / "java_reference.json"
SEGMENTATION_SOURCE_SHA256 = "739a770b866dcc76b8198d795a06f2ae90eecaa437bed50c2fe151bf72bf1fb1"
POWER_TABLE_SOURCE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"


def _load_image(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def test_optimized_tissue_data_terms_match_java_formula_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SOURCE_SHA256
    assert reference["power_table_source_sha256"] == POWER_TABLE_SOURCE_SHA256
    assert reference["shape"] == [7, 6, 5, 3]

    images = _load_image("images")
    field = _load_image("field")
    bias = _load_image("psf_bias")
    x, y, z, channel = np.indices(images.shape, dtype=np.int32)
    expected_images = (((11 * x + 7 * y + 5 * z + 3 * channel) % 101) / 100.0).astype(
        np.float32
    )
    field_mod = ((2 * x + 3 * y + 5 * z + channel) % 16).astype(np.float32)
    expected_field = np.asarray(np.float32(0.85) + field_mod / np.float32(100.0), dtype=np.float32)
    bias_mod = ((5 * x + 3 * y + z + 2 * channel) % 17).astype(np.float32)
    expected_bias = np.asarray(bias_mod / np.float32(1000.0), dtype=np.float32)
    np.testing.assert_array_equal(images, expected_images)
    np.testing.assert_array_equal(field, expected_field)
    np.testing.assert_array_equal(bias, expected_bias)

    parameters = {
        "centroids": np.asarray(reference["centroids"], dtype=np.float32),
        "exponent_by_channel_class": np.asarray(reference["exponents"], dtype=np.float32),
        "modality_weights": np.asarray(reference["modality_weights"], dtype=np.float32),
        "class_variance_weights": np.asarray(
            reference["class_variance_weights"], dtype=np.float32
        ),
        "image_variances": np.asarray(reference["image_variances"], dtype=np.float32),
        "psf_bias": bias,
    }
    without_field = toads_optimized_tissue_data_term(images, **parameters)
    with_field = toads_optimized_tissue_data_term(images, intensity_field=field, **parameters)
    java_without = np.asarray(reference["without_field"], dtype=np.float32).reshape(without_field.shape)
    java_with = np.asarray(reference["with_field"], dtype=np.float32).reshape(with_field.shape)

    np.testing.assert_allclose(without_field, java_without, rtol=0, atol=1e-7)
    np.testing.assert_allclose(with_field, java_with, rtol=0, atol=1e-7)
    assert not np.array_equal(without_field, with_field)
