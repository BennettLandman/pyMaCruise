from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import ImageGeometryError, load_nifti_ras, save_nifti_ras


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
