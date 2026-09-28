from __future__ import annotations

import json
from importlib.metadata import PackageNotFoundError, version

import nibabel as nib
import numpy as np
import pytest

from pymacruise.contracts import WorkflowProfile, validate_workflow_inputs
from pymacruise.workflow import (
    WorkflowConfig,
    create_initial_manifest,
    write_run_manifest,
)


def test_config_rejects_invalid_or_unsupported_values(tmp_path) -> None:
    with pytest.raises(ValueError, match="requires flair_path"):
        WorkflowConfig("plus", "t1.nii.gz", "labels.nii.gz", "mask.nii.gz", tmp_path)
    with pytest.raises(ValueError, match="brain_mask_path must use compressed NIfTI"):
        WorkflowConfig("base", "t1.nii.gz", "labels.nii.gz", "mask.nii", tmp_path)
    with pytest.raises(ValueError, match="Only the source-selected N4"):
        WorkflowConfig(
            "base", "t1.nii.gz", "labels.nii.gz", "mask.nii.gz", tmp_path, bias_correction="N3"
        )


def test_noop_workflow_manifest_captures_inputs_config_and_stage_state(tmp_path) -> None:
    affine = np.diag([1.0, 1.5, 2.0, 1.0])
    t1 = nib.Nifti1Image(np.ones((5, 6, 7), dtype=np.float32), affine)
    labels = np.zeros((5, 6, 7), dtype=np.uint8)
    labels[1:4, 2:5, 2:5] = 4
    maseg = nib.Nifti1Image(labels, affine)
    mask_data = np.zeros((5, 6, 7), dtype=np.uint8)
    mask_data[1:4, 1:5, 1:6] = 255
    brain_mask = nib.Nifti1Image(mask_data, affine)
    t1_path = tmp_path / "t1.nii.gz"
    maseg_path = tmp_path / "maseg.nii.gz"
    mask_path = tmp_path / "brain-mask.nii.gz"
    nib.save(t1, t1_path)
    nib.save(maseg, maseg_path)
    nib.save(brain_mask, mask_path)

    config = WorkflowConfig("base", t1_path, maseg_path, mask_path, tmp_path / "out")
    validated = validate_workflow_inputs(
        WorkflowProfile.BASE, t1, maseg, brain_mask=brain_mask
    )
    manifest = create_initial_manifest(config, validated, source_revision="synthetic-test")
    path = write_run_manifest(config.output_dir / "run.json", manifest)
    saved = json.loads(path.read_text())

    assert saved["workflow_profile"] == "base"
    assert saved["source_revision"] == "synthetic-test"
    assert saved["stages"] == {"input_validation": "passed", "pipeline": "not_started"}
    assert saved["inputs"]["t1"]["output_axis_codes"] == ["R", "A", "S"]
    assert saved["inputs"]["brain_mask"]["output_axis_codes"] == ["R", "A", "S"]
    assert len(saved["inputs"]["t1"]["sha256"]) == 64
    try:
        antspyx_version = version("antspyx")
    except PackageNotFoundError:
        antspyx_version = "not-installed"
    assert saved["dependencies"]["antspyx"] == antspyx_version
    assert len(saved["run_signature"]) == 64


def test_plus_manifest_records_flair_and_signature_changes_with_input(tmp_path) -> None:
    affine = np.eye(4)
    t1 = nib.Nifti1Image(np.ones((4, 4, 4), dtype=np.float32), affine)
    maseg = nib.Nifti1Image(np.zeros((4, 4, 4), dtype=np.uint8), affine)
    flair = nib.Nifti1Image(np.full((5, 5, 5), 2, dtype=np.float32), np.diag([2, 2, 2, 1]))
    brain_mask = nib.Nifti1Image(np.ones((4, 4, 4), dtype=np.uint8), t1.affine)
    paths = [tmp_path / f"{name}.nii.gz" for name in ("t1", "maseg", "brain-mask", "flair")]
    for image, path in zip((t1, maseg, brain_mask, flair), paths, strict=True):
        nib.save(image, path)
    config = WorkflowConfig("plus", paths[0], paths[1], paths[2], tmp_path / "out", paths[3])
    validated = validate_workflow_inputs(
        "plus", t1, maseg, flair, brain_mask=brain_mask
    )
    first = create_initial_manifest(config, validated)
    assert set(first["inputs"]) == {"t1", "maseg", "brain_mask", "flair"}
    changed = nib.Nifti1Image(np.full((4, 4, 4), 3, dtype=np.float32), affine)
    nib.save(changed, paths[0])
    second = create_initial_manifest(config, validated)
    assert second["run_signature"] != first["run_signature"]
