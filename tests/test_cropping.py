from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cropping import (
    CropParameters,
    crop_from_mask,
    crop_image_with_reference_bounds,
    uncrop_image,
)


def test_asymmetric_threshold_bounds_landmark_affine_and_full_uncrop() -> None:
    shape = (17, 19, 21)
    affine = np.diag([0.8, 1.2, 2.0, 1.0])
    affine[:3, 3] = [-30.0, 12.0, -8.0]
    mask_data = np.zeros(shape, dtype=np.uint8)
    mask_data[4:8, 7:12, 9:14] = 1
    # Equal to threshold is excluded by the source's strict comparison.
    mask_data[2, 5, 6] = 0
    target_data = np.fromfunction(
        lambda i, j, k: i + 10 * j + 100 * k, shape, dtype=int
    ).astype(np.float32)
    mask = nib.Nifti1Image(mask_data, affine)
    target = nib.Nifti1Image(target_data, affine)

    cropped_mask, cropped_target, parameters = crop_from_mask(
        mask, target, threshold=0, edge_width=2, background=-5
    )

    assert parameters.bounds == ((2, 9), (5, 13), (7, 15))
    assert parameters.cropped_shape == (8, 9, 9)
    assert cropped_target is not None
    np.testing.assert_allclose(cropped_target.affine, affine @ _translation((2, 5, 7)))
    assert nib.aff2axcodes(cropped_target.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        np.asanyarray(cropped_mask.dataobj), mask_data[2:10, 5:14, 7:16]
    )
    assert np.asanyarray(cropped_target.dataobj)[3, 4, 5] == target_data[5, 9, 12]

    restored = uncrop_image(cropped_target, parameters, legacy_exclusive_max=False)
    expected_restored = np.full(shape, -5, dtype=np.float32)
    expected_restored[2:10, 5:14, 7:16] = target_data[2:10, 5:14, 7:16]
    np.testing.assert_array_equal(np.asanyarray(restored.dataobj), expected_restored)
    np.testing.assert_allclose(restored.affine, affine)

    legacy = np.asanyarray(uncrop_image(cropped_target, parameters).dataobj)
    np.testing.assert_array_equal(legacy[2:9, 5:13, 7:15], target_data[2:9, 5:13, 7:15])
    assert legacy[9, 8, 14] == -5  # Inclusive high corner omitted by JIST.


def test_crop_expands_past_image_with_background_and_jist_uncrop_quirk() -> None:
    shape = (9, 10, 11)
    affine = np.diag([1.0, 1.5, 2.0, 1.0])
    mask_data = np.zeros(shape, dtype=np.uint8)
    mask_data[0:2, 8:10, 5:7] = 1
    target_data = np.arange(np.prod(shape), dtype=np.int16).reshape(shape)
    mask = nib.Nifti1Image(mask_data, affine)
    target = nib.Nifti1Image(target_data, affine)

    _, cropped, parameters = crop_from_mask(mask, target, edge_width=2, background=-9)

    assert cropped is not None
    assert parameters.bounds == ((-2, 3), (6, 11), (3, 8))
    assert cropped.shape == (6, 6, 6)
    values = np.asanyarray(cropped.dataobj)
    assert values[0, 0, 0] == -9
    assert values[2, 2, 2] == target_data[0, 8, 5]
    assert values[-1, -1, -1] == -9
    assert cropped.affine[0, 3] == affine[0, 3] - 2

    legacy = uncrop_image(cropped, parameters)
    restored = np.asanyarray(legacy.dataobj)
    np.testing.assert_array_equal(restored[0:3, 6:9, 3:8], target_data[0:3, 6:9, 3:8])
    assert restored[0, 6, 8] == -9
    assert restored[0, 6, 9] == -9
    np.testing.assert_allclose(legacy.affine, affine)


def test_crop_rejects_empty_masks_and_misaligned_companions() -> None:
    empty = nib.Nifti1Image(np.zeros((5, 6, 7), dtype=np.uint8), np.eye(4))
    with pytest.raises(ValueError, match="empty threshold mask"):
        crop_from_mask(empty)

    mask_data = np.zeros((5, 6, 7), dtype=np.uint8)
    mask_data[2, 2, 2] = 1
    mask = nib.Nifti1Image(mask_data, np.eye(4))
    shifted = np.eye(4)
    shifted[0, 3] = 1
    target = nib.Nifti1Image(np.ones_like(mask_data), shifted)
    with pytest.raises(ValueError, match="same physical grid"):
        crop_from_mask(mask, target)


def test_crop_with_full_reference_bounds_pads_legacy_smaller_extents() -> None:
    affine = np.diag([1.0, 1.5, 2.0, 1.0])
    affine[:3, 3] = (-20.0, 4.0, -10.0)
    values = np.fromfunction(lambda i, j, k: i + 10 * j + 100 * k, (4, 5, 6)).astype(np.float32)
    source = nib.Nifti1Image(values, affine)
    parameters = CropParameters(
        original_shape=(8, 9, 10),
        bounds=((2, 6), (1, 5), (4, 8)),
        background=0.0,
        original_affine=affine,
    )

    result = crop_image_with_reference_bounds(source, parameters)
    expected = np.zeros((5, 5, 5), dtype=np.float32)
    expected[:2, :4, :2] = values[2:4, 1:5, 4:6]

    assert result.shape == (5, 5, 5)
    np.testing.assert_array_equal(np.asanyarray(result.dataobj), expected)
    np.testing.assert_array_equal(result.affine, affine @ _translation((2, 1, 4)))


def test_crop_warns_and_reorients_non_ras_mask_and_target() -> None:
    data = np.zeros((7, 8, 9), dtype=np.uint8)
    data[2:4, 3:5, 4:6] = 1
    non_ras_affine = np.diag([-1.0, -1.0, 1.0, 1.0])
    mask = nib.Nifti1Image(data, non_ras_affine)
    target = nib.Nifti1Image(data.astype(np.float32), non_ras_affine)

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as caught:
        cropped_mask, cropped_target, _ = crop_from_mask(mask, target, edge_width=1)

    assert len(caught) == 2
    assert nib.aff2axcodes(cropped_mask.affine) == ("R", "A", "S")
    assert cropped_target is not None
    assert nib.aff2axcodes(cropped_target.affine) == ("R", "A", "S")


def _translation(offset: tuple[int, int, int]) -> np.ndarray:
    transform = np.eye(4)
    transform[:3, 3] = offset
    return transform
