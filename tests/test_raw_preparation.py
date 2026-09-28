from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.raw_preparation import prepare_raw


def test_prepare_raw_compresses_and_reorients_inputs_to_ras(tmp_path) -> None:
    shape = (7, 8, 9)
    affine = np.diag([1.0, 1.5, 2.0, 1.0])
    t1_data = np.arange(np.prod(shape), dtype=np.int16).reshape(shape)
    labels_data = np.zeros(shape, dtype=np.uint8)
    labels_data[2:5, 2:6, 3:7] = 4
    mask_data = np.zeros(shape, dtype=np.uint8)
    mask_data[1:-1, 1:-1, 1:-1] = 255
    flip = np.eye(4)
    flip[:2, :2] = -np.eye(2)
    flip[0, 3] = shape[0] - 1
    flip[1, 3] = shape[1] - 1
    t1_path = tmp_path / "source-t1.nii"
    labels_path = tmp_path / "source-maseg.nii.gz"
    mask_path = tmp_path / "source-mask.nii.gz"
    nib.save(nib.Nifti1Image(t1_data[::-1, ::-1, :], affine @ flip), t1_path)
    nib.save(nib.Nifti1Image(labels_data[::-1, ::-1, :], affine @ flip), labels_path)
    nib.save(nib.Nifti1Image(mask_data[::-1, ::-1, :], affine @ flip), mask_path)

    with pytest.warns(UserWarning) as warnings_seen:
        result = prepare_raw(t1_path, labels_path, mask_path, tmp_path / "prepared")
    assert len(warnings_seen) == 3
    assert result.input_orientations == {
        "t1": ("L", "P", "S"),
        "maseg": ("L", "P", "S"),
        "brain_mask": ("L", "P", "S"),
    }
    assert result.t1_path.name == "T1.nii.gz"
    assert result.maseg_path.name == "Seg.nii.gz"
    assert result.brain_mask_path.name == "BrainMask.nii.gz"
    t1 = nib.load(result.t1_path)
    maseg = nib.load(result.maseg_path)
    mask = nib.load(result.brain_mask_path)
    assert nib.aff2axcodes(t1.affine) == nib.aff2axcodes(maseg.affine) == nib.aff2axcodes(mask.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asanyarray(t1.dataobj), t1_data)
    np.testing.assert_array_equal(np.asanyarray(maseg.dataobj), labels_data)
    np.testing.assert_array_equal(np.asanyarray(mask.dataobj), mask_data)
    np.testing.assert_allclose(t1.affine, affine, atol=1e-5)
    np.testing.assert_allclose(maseg.affine, affine, atol=1e-5)
    np.testing.assert_allclose(mask.affine, affine, atol=1e-5)


def test_prepare_raw_rejects_stale_outputs_and_non_nifti_inputs(tmp_path) -> None:
    data = np.ones((4, 4, 4), dtype=np.float32)
    affine = np.eye(4)
    t1 = tmp_path / "t1.nii.gz"
    labels = tmp_path / "labels.nii.gz"
    mask = tmp_path / "mask.nii.gz"
    nib.save(nib.Nifti1Image(data, affine), t1)
    nib.save(nib.Nifti1Image(np.zeros((4, 4, 4), dtype=np.uint8), affine), labels)
    mask_data = np.zeros((4, 4, 4), dtype=np.uint8)
    mask_data[1:3, 1:3, 1:3] = 1
    nib.save(nib.Nifti1Image(mask_data, affine), mask)
    output_dir = tmp_path / "prepared"
    prepare_raw(t1, labels, mask, output_dir)
    with pytest.raises(FileExistsError, match="overwrite=True"):
        prepare_raw(t1, labels, mask, output_dir)
    with pytest.raises(ValueError, match="convert other formats first"):
        prepare_raw(tmp_path / "image.dcm", labels, mask, tmp_path / "other")
