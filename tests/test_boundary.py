from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.boundary import (
    clean_csf_boundary,
    clean_csf_boundary_image,
    merge_csf_boundary,
    merge_csf_boundary_images,
    multi_atlas_csf_boundary,
    multi_atlas_csf_boundary_image,
)


def test_default_boundary_is_the_exact_two_label_interface() -> None:
    segmentation = np.zeros((5, 5, 5), dtype=np.int16)
    segmentation[:2] = 100
    segmentation[2:] = 101

    expected = np.zeros_like(segmentation, dtype=np.uint8)
    expected[1:3, 1:4, 1:4] = 1
    np.testing.assert_array_equal(multi_atlas_csf_boundary(segmentation), expected)


def test_overlooked_labels_are_removed_and_low_label_contacts_suppress_boundary() -> None:
    segmentation = np.zeros((7, 7, 7), dtype=np.int16)
    segmentation[:3] = 100
    segmentation[3:] = 101
    segmentation[3, 2, 2] = 120
    segmentation[3, 5, 5] = 5

    boundary = multi_atlas_csf_boundary(segmentation)
    assert boundary[3, 2, 2] == 0  # overlooked 120 was zeroed
    assert boundary[3, 4, 5] == 0  # a low label suppresses this interface voxel
    assert boundary[3, 3, 3] == 1


def test_boundary_handles_empty_and_small_volumes() -> None:
    assert not multi_atlas_csf_boundary(np.zeros((4, 4, 4), dtype=np.uint8)).any()
    assert not multi_atlas_csf_boundary(np.ones((2, 4, 4), dtype=np.uint8) * 100).any()
    with pytest.raises(ValueError, match="3D"):
        multi_atlas_csf_boundary(np.ones((4, 4), dtype=np.uint8))


def test_merge_modes_match_source_operations_and_reject_bad_shapes() -> None:
    boundary = np.array([[[0, 1], [1, 0]]], dtype=np.uint8)
    csf = np.array([[[0.6, 0.4], [1.1, -0.2]]], dtype=np.float32)
    merged = merge_csf_boundary(boundary, csf, mode="boundary_as_one")
    np.testing.assert_allclose(merged, [[[0.6, 1.0], [1.0, -0.2]]])

    speed = merge_csf_boundary(boundary, csf, mode="speed_function")
    np.testing.assert_allclose(speed, [[[0.6, 0.1], [0.1, -0.2]]])
    with pytest.raises(ValueError, match="shapes differ"):
        merge_csf_boundary(boundary, csf[..., :1], mode="boundary_as_one")
    with pytest.raises(ValueError, match="Unsupported"):
        merge_csf_boundary(boundary, csf, mode="unknown")  # type: ignore[arg-type]


def test_nifti_boundary_and_merge_keep_ras_geometry() -> None:
    affine = np.diag([1.5, 2.0, 2.5, 1.0])
    labels = np.zeros((5, 5, 5), dtype=np.int16)
    labels[:2] = 100
    labels[2:] = 101
    segmentation = nib.Nifti1Image(labels, affine)
    boundary = multi_atlas_csf_boundary_image(segmentation)
    assert nib.aff2axcodes(boundary.affine) == ("R", "A", "S")
    assert boundary.get_data_dtype() == np.dtype(np.uint8)
    np.testing.assert_array_equal(boundary.affine, affine)

    csf = nib.Nifti1Image(np.full(labels.shape, 0.25, dtype=np.float32), affine)
    merged = merge_csf_boundary_images(boundary, csf, mode="boundary_as_one")
    assert nib.aff2axcodes(merged.affine) == ("R", "A", "S")
    assert merged.get_data_dtype() == np.dtype(np.float32)
    assert np.all(np.asanyarray(merged.dataobj)[np.asanyarray(boundary.dataobj) > 0] == 1.0)


def test_image_merge_rejects_mismatched_physical_grids() -> None:
    first = nib.Nifti1Image(np.zeros((3, 3, 3), dtype=np.uint8), np.eye(4))
    shifted = np.eye(4)
    shifted[0, 3] = 1
    second = nib.Nifti1Image(np.zeros((3, 3, 3), dtype=np.float32), shifted)
    with pytest.raises(ValueError, match="same physical grid"):
        merge_csf_boundary_images(first, second, mode="boundary_as_one")


@pytest.mark.parametrize(
    ("mode", "threshold"),
    [("toads_wm", 0.01), ("updated_wm", 25.5)],
)
def test_cleanup_uses_strict_source_threshold(mode: str, threshold: float) -> None:
    boundary = np.ones((1, 1, 4), dtype=np.uint8)
    wm = np.array([[[threshold, np.nextafter(threshold, np.inf), 0.0, threshold]]])
    cleaned = clean_csf_boundary(boundary, wm, mode=mode)  # type: ignore[arg-type]
    np.testing.assert_array_equal(cleaned, [[[1, 0, 1, 1]]])


def test_cleanup_nifti_retains_ras_affine() -> None:
    affine = np.diag([1.0, 1.25, 2.0, 1.0])
    boundary = nib.Nifti1Image(np.ones((3, 3, 3), dtype=np.uint8), affine)
    wm_data = np.zeros((3, 3, 3), dtype=np.float32)
    wm_data[1, 1, 1] = 0.02
    wm = nib.Nifti1Image(wm_data, affine)
    cleaned = clean_csf_boundary_image(boundary, wm, mode="toads_wm")
    assert nib.aff2axcodes(cleaned.affine) == ("R", "A", "S")
    assert np.asanyarray(cleaned.dataobj)[1, 1, 1] == 0
