from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_psf_ridge_bias


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_psf_bias"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
RIDGE = np.float32(0.37)


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def _expected_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = (9, 8, 7, 2)
    x, y, z, channel = np.indices(shape, dtype=np.int32)
    images = (((13 * x + 11 * y + 7 * z + 5 * channel) % 113) / 112.0).astype(np.float32)
    field_remainder = ((3 * x + 5 * y + 2 * z + channel) % 21).astype(np.float32)
    field = np.asarray(np.float32(0.8) + field_remainder / np.float32(100.0), dtype=np.float32)
    segmentation = np.broadcast_to((x[..., 0] >= 4).astype(np.int16), shape[:3]).copy()
    segmentation[2, 3, 3] = 2
    segmentation[6, 1:4, 2:5] = 0
    return images, field, segmentation


def test_psf_ridge_bias_matches_java_source_formula() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["shape"] == [9, 8, 7, 2]
    assert reference["ridge"] == 0.37

    images, field, segmentation = (_load("images"), _load("field"), _load("segmentation"))
    expected_images, expected_field, expected_segmentation = _expected_inputs()
    np.testing.assert_array_equal(images, expected_images)
    np.testing.assert_array_equal(field, expected_field)
    np.testing.assert_array_equal(segmentation, expected_segmentation)

    without_field = toads_psf_ridge_bias(images, segmentation, RIDGE)
    with_field = toads_psf_ridge_bias(images, segmentation, RIDGE, intensity_field=field)
    java_without = np.asarray(reference["without_field"], dtype=np.float32).reshape(without_field.shape)
    java_with = np.asarray(reference["with_field"], dtype=np.float32).reshape(with_field.shape)
    np.testing.assert_array_equal(without_field, java_without)
    np.testing.assert_array_equal(with_field, java_with)

    assert np.count_nonzero(without_field) > 0
    assert not np.array_equal(without_field, with_field)
    assert np.count_nonzero(without_field[0]) == 0
    assert np.count_nonzero(without_field[:, -1]) == 0
    assert np.count_nonzero(without_field[:, :, 0]) == 0


def test_psf_ridge_bias_zero_ridge_and_nonboundary_return_zero() -> None:
    images = np.ones((5, 5, 5, 1), dtype=np.float32)
    uniform = np.zeros((5, 5, 5), dtype=np.uint8)
    assert not np.any(toads_psf_ridge_bias(images, uniform, 0.7))
    segmented = uniform.copy()
    segmented[3:, :, :] = 1
    assert not np.any(toads_psf_ridge_bias(images, segmented, 0.0))
    bias = toads_psf_ridge_bias(images, segmented, 0.7)
    assert np.all(bias[1, 1, 1] == 0.0)
    assert np.any(bias[2, 2, 2] != 0.0)


@pytest.mark.parametrize(
    ("images", "segmentation", "ridge", "field"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3), dtype=np.uint8), 1.0, None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((2, 3, 3), dtype=np.uint8), 1.0, None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3), dtype=np.float32), 1.0, None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3), dtype=np.uint8), np.nan, None),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3), dtype=np.uint8), 1.0, np.zeros((3, 3, 3))),
    ],
)
def test_psf_ridge_bias_rejects_invalid_inputs(images, segmentation, ridge, field) -> None:
    with pytest.raises(ValueError):
        toads_psf_ridge_bias(images, segmentation, ridge, intensity_field=field)
