from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.contracts import (
    WorkflowInputError,
    WorkflowProfile,
    validate_workflow_files,
    validate_workflow_inputs,
)


def _aligned_inputs() -> tuple[nib.Nifti1Image, nib.Nifti1Image]:
    affine = np.diag([1.0, 1.2, 1.5, 1.0])
    t1 = nib.Nifti1Image(np.ones((6, 7, 8), dtype=np.float32), affine)
    labels = np.zeros((6, 7, 8), dtype=np.uint8)
    labels[1:3, 2:5, 3:7] = 4
    labels[3:5, 2:5, 3:7] = 100
    maseg = nib.Nifti1Image(labels, affine)
    return t1, maseg


def test_base_contract_reorients_images_and_records_original_axes() -> None:
    t1, maseg = _aligned_inputs()
    flip = np.diag([-1.0, -1.0, 1.0, 1.0])
    flip[0, 3] = t1.shape[0] - 1
    flip[1, 3] = t1.shape[1] - 1
    non_ras_labels = nib.Nifti1Image(np.asanyarray(maseg.dataobj)[::-1, ::-1, :], t1.affine @ flip)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = validate_workflow_inputs("base", t1, non_ras_labels)
    assert result.profile is WorkflowProfile.BASE
    assert result.input_orientations["maseg"] == ("L", "P", "S")
    assert nib.aff2axcodes(result.maseg.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(result.maseg.get_fdata(), maseg.get_fdata())


def test_plus_requires_flair_and_accepts_a_distinct_input_grid() -> None:
    t1, maseg = _aligned_inputs()
    with pytest.raises(WorkflowInputError, match="requires a FLAIR"):
        validate_workflow_inputs("plus", t1, maseg)
    with pytest.raises(WorkflowInputError, match="does not accept"):
        validate_workflow_inputs("base", t1, maseg, t1)

    flair_affine = np.diag([0.8, 0.8, 2.0, 1.0])
    flair = nib.Nifti1Image(np.ones((9, 10, 11), dtype=np.float32), flair_affine)
    result = validate_workflow_inputs("plus", t1, maseg, flair)
    assert result.flair is flair
    assert result.profile is WorkflowProfile.PLUS


@pytest.mark.parametrize("problem", ["shifted", "non_integer", "infinite", "4d"])
def test_invalid_maseg_and_geometry_fail_with_actionable_errors(problem: str) -> None:
    t1, maseg = _aligned_inputs()
    if problem == "shifted":
        affine = maseg.affine.copy()
        affine[0, 3] = 2
        maseg = nib.Nifti1Image(np.asanyarray(maseg.dataobj), affine)
        message = "not aligned"
    elif problem == "non_integer":
        data = np.asanyarray(maseg.dataobj, dtype=np.float32)
        data[0, 0, 0] = 0.5
        maseg = nib.Nifti1Image(data, maseg.affine)
        message = "integer-valued"
    elif problem == "infinite":
        data = np.asanyarray(maseg.dataobj, dtype=np.float32)
        data[0, 0, 0] = np.inf
        maseg = nib.Nifti1Image(data, maseg.affine)
        message = "NaN or infinite"
    else:
        maseg = nib.Nifti1Image(np.zeros((6, 7, 8, 1)), maseg.affine)
        message = "3D single-volume"
    with pytest.raises(WorkflowInputError, match=message):
        validate_workflow_inputs("base", t1, maseg)


def test_file_contract_requires_compressed_nifti_and_reports_missing_path(tmp_path) -> None:
    with pytest.raises(WorkflowInputError, match="Missing T1 input"):
        validate_workflow_files("base", tmp_path / "missing.nii.gz", tmp_path / "labels.nii.gz")

    t1, maseg = _aligned_inputs()
    t1_path = tmp_path / "t1.nii.gz"
    maseg_path = tmp_path / "maseg.nii.gz"
    nib.save(t1, t1_path)
    nib.save(maseg, maseg_path)
    assert validate_workflow_files("base", t1_path, maseg_path).profile is WorkflowProfile.BASE
    wrong_suffix = tmp_path / "t1.nii"
    wrong_suffix.write_bytes(b"placeholder")
    with pytest.raises(WorkflowInputError, match=".nii.gz"):
        validate_workflow_files("base", wrong_suffix, maseg_path)
