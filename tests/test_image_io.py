from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import (
    ImageGeometryError,
    is_axis_aligned_ras_grid,
    is_orthogonal_ras_grid,
    load_nifti_ras,
    regrid_to_orthogonal_ras,
    save_nifti_ras,
)


def test_ras_nifti_round_trip_preserves_values_geometry_and_header(tmp_path) -> None:
    data = np.arange(5 * 6 * 7, dtype=np.int16).reshape((5, 6, 7))
    affine = np.array(
        [
            [0.0, -2.0, 0.0, 10.0],
            [1.5, 0.0, 0.0, -4.0],
            [0.0, 0.0, 3.0, 8.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    image = nib.Nifti1Image(data, affine)
    image.header.set_xyzt_units("mm", "sec")
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    input_path = tmp_path / "input.nii.gz"
    output_path = tmp_path / "nested" / "output.nii.gz"
    nib.save(image, input_path)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        loaded = load_nifti_ras(input_path)
    assert nib.aff2axcodes(loaded.affine) == ("R", "A", "S")
    for source_index in ((0, 0, 0), (4, 5, 6), (2, 3, 1)):
        world = nib.affines.apply_affine(image.affine, source_index)
        ras_index = np.rint(nib.affines.apply_affine(np.linalg.inv(loaded.affine), world)).astype(int)
        assert np.asarray(loaded.dataobj)[tuple(ras_index)] == np.asarray(image.dataobj)[source_index]
        np.testing.assert_allclose(
            nib.affines.apply_affine(loaded.affine, ras_index), world, atol=1e-5
        )
    save_nifti_ras(loaded, output_path)
    written = nib.load(output_path)
    assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
    assert written.get_data_dtype() == np.dtype(np.int16)
    assert written.header.get_xyzt_units() == ("mm", "sec")
    np.testing.assert_array_equal(np.asanyarray(written.dataobj), np.asanyarray(loaded.dataobj))
    np.testing.assert_allclose(written.affine, loaded.affine, atol=1e-5)
    assert written.header["sform_code"] == 1


def test_load_rejects_conflicting_forms_and_non_nifti_suffixes(tmp_path) -> None:
    image = nib.Nifti1Image(np.ones((3, 4, 5), dtype=np.float32), np.eye(4))
    different_qform = np.eye(4)
    different_qform[0, 3] = 5
    image.set_qform(different_qform, code=1)
    image.set_sform(np.eye(4), code=1)
    path = tmp_path / "ambiguous.nii.gz"
    nib.save(image, path)
    with pytest.raises(ImageGeometryError, match="conflicting"):
        load_nifti_ras(path)
    with pytest.raises(ValueError, match=".nii.gz"):
        load_nifti_ras(tmp_path / "image.nii")

    malformed = tmp_path / "malformed.nii.gz"
    malformed.write_bytes(b"not a NIfTI file")
    with pytest.raises(ValueError, match="Cannot read NIfTI"):
        load_nifti_ras(malformed)


def test_scaled_voxel_values_survive_decode_and_resave(tmp_path) -> None:
    raw = np.arange(24, dtype=np.int16).reshape((2, 3, 4))
    image = nib.Nifti1Image(raw, np.eye(4))
    image.header.set_slope_inter(2.0, -3.0)
    source = tmp_path / "scaled.nii.gz"
    destination = tmp_path / "decoded.nii.gz"
    nib.save(image, source)
    loaded = load_nifti_ras(source)
    expected = np.asanyarray(loaded.dataobj).copy()
    save_nifti_ras(loaded, destination)
    resaved = nib.load(destination)
    np.testing.assert_array_equal(np.asanyarray(resaved.dataobj), expected)


def test_oblique_ras_inputs_regrid_to_one_orthogonal_grid_with_role_interpolation() -> None:
    angle = np.deg2rad(12.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    affine = np.eye(4)
    affine[:3, :3] = rotation @ np.diag([1.0, 1.25, 1.5])
    affine[:3, 3] = [-8.0, -7.0, -6.0]
    shape = (12, 11, 10)
    x, y, z = np.indices(shape, dtype=np.float32)
    t1 = nib.Nifti1Image((x + 2 * y + 3 * z).astype(np.float32), affine)
    mask_values = np.zeros(shape, dtype=np.uint8)
    mask_values[2:9, 2:9, 2:8] = 1
    mask = nib.Nifti1Image(mask_values, affine)
    label_values = np.zeros(shape, dtype=np.int16)
    label_values[3:8, 3:8, 3:7] = 100
    maseg = nib.Nifti1Image(label_values, affine)

    target_t1 = regrid_to_orthogonal_ras(t1, interpolation_order=1)
    target_mask = regrid_to_orthogonal_ras(mask, reference=target_t1, interpolation_order=0)
    target_maseg = regrid_to_orthogonal_ras(maseg, reference=target_t1, interpolation_order=0)

    assert is_axis_aligned_ras_grid(target_t1)
    assert all(image.shape == target_t1.shape for image in (target_mask, target_maseg))
    assert all(np.allclose(image.affine, target_t1.affine) for image in (target_mask, target_maseg))
    np.testing.assert_array_equal(np.unique(np.asarray(target_mask.dataobj)), [0, 1])
    np.testing.assert_array_equal(np.unique(np.asarray(target_maseg.dataobj)), [0, 100])
    assert np.isfinite(np.asarray(target_t1.dataobj)).all()


def test_orthogonal_ras_grid_recognizes_oblique_without_regridding() -> None:
    angle = np.deg2rad(17.0)
    affine = np.eye(4)
    affine[:3, :3] = np.asarray(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    ) @ np.diag((1.0, 1.2, 1.5))
    image = nib.Nifti1Image(np.zeros((4, 5, 6), dtype=np.float32), affine)

    assert is_orthogonal_ras_grid(image)
    assert not is_axis_aligned_ras_grid(image)


def test_axis_aligned_ras_regrid_is_a_noop_and_rejects_other_interpolators() -> None:
    image = nib.Nifti1Image(np.arange(24, dtype=np.int16).reshape((2, 3, 4)), np.eye(4))
    assert is_axis_aligned_ras_grid(image)
    assert regrid_to_orthogonal_ras(image) is image
    with pytest.raises(ValueError, match="interpolation order"):
        regrid_to_orthogonal_ras(image, interpolation_order=2)
