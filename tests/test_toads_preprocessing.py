from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_preprocessing import prepare_toads_intensity, robust_intensity_bounds


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
