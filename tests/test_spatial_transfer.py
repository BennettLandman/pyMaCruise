from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.spatial_transfer import transfer_to_reference_space


def test_known_translation_transfers_labels_by_physical_nearest_neighbor() -> None:
    data = np.zeros((5, 5, 5), dtype=np.uint8)
    data[2, 2, 2] = 7
    source = nib.Nifti1Image(data, np.eye(4))
    translated_reference = np.eye(4)
    translated_reference[0, 3] = 1.0
    reference = nib.Nifti1Image(np.zeros_like(data), translated_reference)

    output = transfer_to_reference_space(source, reference)
    expected = np.zeros_like(data)
    expected[1, 2, 2] = 7
    np.testing.assert_array_equal(np.asanyarray(output.dataobj), expected)
    np.testing.assert_array_equal(output.affine, translated_reference)
    assert output.get_data_dtype() == np.dtype(np.uint8)


def test_known_rotation_transfers_continuous_values_without_linear_mixing() -> None:
    shape = (5, 5, 5)
    source_data = np.fromfunction(
        lambda i, j, k: i + 10 * j + 100 * k, shape, dtype=int
    ).astype(np.float32)
    source_affine = np.diag([1.0, 1.0, 1.0, 1.0])
    source = nib.Nifti1Image(source_data, source_affine)

    rotated = np.eye(4)
    rotated[:3, :3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    rotated[0, 3] = shape[0] - 1
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), rotated)
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        output = transfer_to_reference_space(source, reference)

    # The reference's 90-degree rotation about the volume center maps output
    # index (i,j,k) to source index (N-1-j, i, k). Compare only this exact
    # integer-lattice geometry; there is no interpolation between intensities.
    target_ras = nib.as_closest_canonical(reference)
    inverse_source = np.linalg.inv(source_affine)
    expected = np.zeros(target_ras.shape, dtype=np.float32)
    for index in np.ndindex(target_ras.shape):
        world = target_ras.affine @ np.array([*index, 1.0])
        source_index = np.rint((inverse_source @ world)[:3]).astype(int)
        if np.all(source_index >= 0) and np.all(source_index < shape):
            expected[index] = source_data[tuple(source_index)]

    np.testing.assert_array_equal(np.asanyarray(output.dataobj), expected)
    assert output.shape == target_ras.shape
    assert nib.aff2axcodes(output.affine) == ("R", "A", "S")
    assert output.get_data_dtype() == np.dtype(np.float32)


def test_non_ras_source_and_reference_warn_and_output_is_ras() -> None:
    data = np.zeros((4, 4, 4), dtype=np.uint8)
    data[1, 1, 1] = 3
    source = nib.Nifti1Image(data, np.diag([-1.0, -1.0, 1.0, 1.0]))
    reference = nib.Nifti1Image(np.zeros_like(data), np.eye(4))
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        output = transfer_to_reference_space(source, reference)
    assert nib.aff2axcodes(output.affine) == ("R", "A", "S")


def test_transfer_rejects_nonfinite_input_and_non_3d() -> None:
    bad = np.ones((3, 3, 3), dtype=np.float32)
    bad[0, 0, 0] = np.nan
    source = nib.Nifti1Image(bad, np.eye(4))
    reference = nib.Nifti1Image(np.zeros_like(bad), np.eye(4))
    with pytest.raises(ValueError, match="NaN"):
        transfer_to_reference_space(source, reference)

    four_d = nib.Nifti1Image(np.ones((3, 3, 3, 2), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="3D"):
        transfer_to_reference_space(four_d, reference)
