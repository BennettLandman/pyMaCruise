from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.volume_ops import VolumeGeometryError, calculate_nifti_files, calculate_volume


def test_voxel_operation_reorients_non_ras_inputs_with_warning() -> None:
    data = np.arange(5 * 6 * 7, dtype=np.int16).reshape((5, 6, 7))
    ras_affine = np.diag([1.0, 2.0, 3.0, 1.0])
    image_a = nib.Nifti1Image(data, ras_affine)

    flip = np.array(
        [
            [-1.0, 0.0, 0.0, data.shape[0] - 1],
            [0.0, -1.0, 0.0, data.shape[1] - 1],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    image_b = nib.Nifti1Image(data[::-1, ::-1, :].copy(), ras_affine @ flip)
    assert nib.aff2axcodes(image_b.affine) == ("L", "P", "S")

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = calculate_volume(image_a, image_b, "ADD")

    assert nib.aff2axcodes(result.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asanyarray(result.dataobj), data.astype(np.float32) * 2)


def test_voxel_operation_rejects_misaligned_ras_images() -> None:
    data = np.ones((4, 5, 6), dtype=np.float32)
    first = nib.Nifti1Image(data, np.eye(4))
    shifted_affine = np.eye(4)
    shifted_affine[0, 3] = 1.0
    second = nib.Nifti1Image(data, shifted_affine)

    with pytest.raises(VolumeGeometryError, match="affines differ"):
        calculate_volume(first, second, "ADD")


def test_nifti_file_operator_enforces_compressed_ras_output(tmp_path) -> None:
    affine = np.diag([-1.0, 2.0, 3.0, 1.0])
    image_a = nib.Nifti1Image(np.ones((4, 5, 6), dtype=np.uint8), affine)
    image_b = nib.Nifti1Image(np.full((4, 5, 6), 2, dtype=np.uint8), affine)
    input_a = tmp_path / "a.nii.gz"
    input_b = tmp_path / "b.nii.gz"
    output = tmp_path / "result.nii.gz"
    nib.save(image_a, input_a)
    nib.save(image_b, input_b)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        calculate_nifti_files(str(input_a), str(input_b), str(output), "ADD")

    result = nib.load(output)
    assert nib.aff2axcodes(result.affine) == ("R", "A", "S")
    assert np.all(np.asanyarray(result.dataobj) == 3)
    with pytest.raises(ValueError, match=".nii.gz"):
        calculate_nifti_files(str(input_a), str(input_b), str(tmp_path / "out.nii"), "ADD")
