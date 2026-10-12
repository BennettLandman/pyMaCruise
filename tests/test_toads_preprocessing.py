from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_preprocessing import (
    prepare_toads_intensity,
    robust_intensity_bounds,
    toads_ras_geometry,
)


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "toads_preprocessing"


def test_toads_robust_normalization_matches_java_reference() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    mask = nib.load(FIXTURES / "brain_mask.nii.gz")
    reference = json.loads((FIXTURES / "java_reference.json").read_text())

    result = prepare_toads_intensity(t1, mask)

    assert nib.aff2axcodes(result.image.affine) == ("R", "A", "S")
    assert result.image.shape == tuple(reference["shape"])
    assert result.robust_minimum == np.float32(reference["low"])
    assert result.robust_maximum == np.float32(reference["high"])
    assert result.masked_voxels == 594
    actual = np.asanyarray(result.image.dataobj).ravel(order="F")
    np.testing.assert_array_equal(actual, np.asarray(reference["normalized"], dtype=np.float32))
    assert np.all(actual.reshape(reference["shape"], order="F")[
        np.asanyarray(mask.dataobj) == 0
    ] == 0)


def test_toads_robust_bounds_reject_invalid_arrays_and_ranges() -> None:
    with pytest.raises(ValueError, match="nonempty 3D"):
        robust_intensity_bounds(np.ones((3, 4), dtype=np.float32))
    with pytest.raises(ValueError, match="NaN or infinite"):
        robust_intensity_bounds(np.full((3, 4, 5), np.nan, dtype=np.float32))
    with pytest.raises(ValueError, match="fraction"):
        robust_intensity_bounds(np.ones((3, 4, 5), dtype=np.float32), fraction=-0.1)


def test_toads_mask_grid_mismatch_is_rejected() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    mask = nib.load(FIXTURES / "brain_mask.nii.gz")
    affine = mask.affine.copy()
    affine[0, 3] += 1
    shifted_mask = nib.Nifti1Image(np.asanyarray(mask.dataobj), affine)
    with pytest.raises(ValueError, match="not aligned"):
        prepare_toads_intensity(t1, shifted_mask)


def test_toads_non_ras_inputs_warn_and_preserve_normalized_values() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    mask = nib.load(FIXTURES / "brain_mask.nii.gz")
    flip = np.eye(4)
    flip[0, 0] = -1
    flip[1, 1] = -1
    flip[0, 3] = t1.shape[0] - 1
    flip[1, 3] = t1.shape[1] - 1
    non_ras_t1 = nib.Nifti1Image(np.asanyarray(t1.dataobj)[::-1, ::-1, :], t1.affine @ flip)
    non_ras_mask = nib.Nifti1Image(
        np.asanyarray(mask.dataobj)[::-1, ::-1, :], mask.affine @ flip
    )
    expected = prepare_toads_intensity(t1, mask)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        actual = prepare_toads_intensity(non_ras_t1, non_ras_mask)

    assert nib.aff2axcodes(actual.image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(actual.image.dataobj, expected.image.dataobj)


def test_toads_ras_geometry_maps_mipav_directions_to_signed_spacing() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    geometry = toads_ras_geometry(t1)

    assert nib.aff2axcodes(geometry.image.affine) == ("R", "A", "S")
    assert geometry.source_orientation_codes == ("R", "A", "S")
    assert geometry.image_orientation == 0  # MIPAV FileInfoBase.AXIAL
    assert geometry.axis_orientation_codes == (2, 3, 5)  # L2R, P2A, I2S
    assert geometry.orientation_signs == (-1.0, -1.0, 1.0)
    zooms = nib.affines.voxel_sizes(t1.affine)
    np.testing.assert_array_equal(
        geometry.signed_image_resolutions_mm,
        np.asarray((-zooms[0], -zooms[1], zooms[2]), dtype=np.float64),
    )
    assert geometry.initial_transform == (0.0,) * 6


def test_toads_ras_geometry_registration_is_world_origin_independent() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    shifted_affine = t1.affine.copy()
    shifted_affine[:3, 3] += np.asarray((31.0, -17.0, 8.5))
    shifted = nib.Nifti1Image(np.asanyarray(t1.dataobj), shifted_affine)

    original = toads_ras_geometry(t1)
    translated = toads_ras_geometry(shifted)

    np.testing.assert_array_equal(
        translated.signed_image_resolutions_mm,
        original.signed_image_resolutions_mm,
    )
    assert translated.initial_transform == original.initial_transform
    # NIfTI world origins remain in output metadata, but the selected source
    # mapping receives only dimensions and signed spacings.
    np.testing.assert_array_equal(
        translated.image.affine[:3, 3], original.image.affine[:3, 3] + (31.0, -17.0, 8.5)
    )


def test_toads_ras_geometry_reorients_before_deriving_mipav_directions() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    flip = np.diag((-1.0, -1.0, 1.0, 1.0))
    flip[0, 3] = t1.shape[0] - 1
    flip[1, 3] = t1.shape[1] - 1
    non_ras = nib.Nifti1Image(
        np.asanyarray(t1.dataobj)[::-1, ::-1, :], t1.affine @ flip
    )

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        geometry = toads_ras_geometry(non_ras)

    assert geometry.source_orientation_codes == ("L", "P", "S")
    assert nib.aff2axcodes(geometry.image.affine) == ("R", "A", "S")
    assert geometry.orientation_signs == (-1.0, -1.0, 1.0)
    np.testing.assert_array_equal(geometry.image.dataobj, t1.dataobj)


@pytest.mark.parametrize("target_codes", [("A", "R", "S"), ("I", "L", "A")])
def test_toads_ras_geometry_reorients_permuted_axis_orders(
    target_codes: tuple[str, str, str],
) -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    transform = nib.orientations.ornt_transform(
        nib.orientations.axcodes2ornt(("R", "A", "S")),
        nib.orientations.axcodes2ornt(target_codes),
    )
    data = nib.orientations.apply_orientation(np.asanyarray(t1.dataobj), transform)
    affine = t1.affine @ nib.orientations.inv_ornt_aff(transform, t1.shape)
    non_ras = nib.Nifti1Image(data, affine)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        geometry = toads_ras_geometry(non_ras)

    assert geometry.source_orientation_codes == target_codes
    assert nib.aff2axcodes(geometry.image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(geometry.image.dataobj, t1.dataobj)
    np.testing.assert_array_equal(
        geometry.signed_image_resolutions_mm,
        toads_ras_geometry(t1).signed_image_resolutions_mm,
    )


def test_toads_ras_geometry_preserves_oblique_orthogonal_affines() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    affine = t1.affine.copy()
    angle = np.deg2rad(11.0)
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    affine[:3, :3] = rotation @ affine[:3, :3]
    oblique = nib.Nifti1Image(np.asanyarray(t1.dataobj), affine)
    geometry = toads_ras_geometry(oblique)

    assert nib.aff2axcodes(geometry.image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        geometry.image.affine[:3, :3], affine[:3, :3], rtol=0, atol=1e-5
    )
    np.testing.assert_array_equal(
        geometry.signed_image_resolutions_mm,
        np.asarray((-1.0, -1.0, 1.0), dtype=np.float64),
    )


def test_toads_ras_geometry_rejects_sheared_affines() -> None:
    t1 = nib.load(FIXTURES / "t1_n4.nii.gz")
    affine = t1.affine.copy()
    affine[0, 1] = 0.1
    sheared = nib.Nifti1Image(np.asanyarray(t1.dataobj), affine)
    with pytest.raises(ValueError, match="orthogonal RAS voxel axes"):
        toads_ras_geometry(sheared)
