from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

import pymacruise.plus_workflow as plus_workflow
from pymacruise.registration import RegistrationConfig


FIXTURES = Path(__file__).parent / "data/synthetic/registration"


def _image(name: str) -> nib.Nifti1Image:
    image = nib.load(str(FIXTURES / f"{name}.nii.gz"))
    assert isinstance(image, nib.Nifti1Image)
    return image


def test_plus_flair_registration_checkpoint_persists_and_resumes(tmp_path: Path, monkeypatch) -> None:
    config = RegistrationConfig(
        transform_type="Rigid",
        metric="mattes",
        sampling_rate=1.0,
        iterations=(60, 40, 20, 0),
        shrink_factors=(6, 4, 2, 1),
        smoothing_sigmas=(3, 2, 1, 0),
    )
    output_dir = tmp_path / "plus-registration"
    first = plus_workflow.run_plus_flair_registration_checkpoint(
        FIXTURES / "fixed_image.nii.gz",
        FIXTURES / "moving_flair.nii.gz",
        output_dir,
        config=config,
    )
    fixed = _image("fixed_image")
    assert not first.resumed
    assert first.warped_flair_path.name == "flair_in_t1.nii.gz"
    assert first.forward_transforms and first.inverse_transforms
    assert all(path.is_file() and path.is_relative_to(output_dir) for path in (
        *first.forward_transforms, *first.inverse_transforms
    ))
    assert first.warped_flair.shape == fixed.shape
    np.testing.assert_allclose(first.warped_flair.affine, fixed.affine, rtol=0, atol=1e-6)
    assert nib.aff2axcodes(first.warped_flair.affine) == ("R", "A", "S")
    assert np.isfinite(np.asarray(first.warped_flair.dataobj)).all()
    manifest = json.loads(first.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["registration_backend"] == "ANTsPy"
    assert manifest["spatial_convention"] == "RAS"
    assert "Aladin numerical equivalence is unvalidated" in manifest["equivalence_note"]

    def unexpected_registration(*args, **kwargs):
        raise AssertionError("valid resume should not repeat ANTs registration")

    monkeypatch.setattr(plus_workflow, "register_flair_to_t1", unexpected_registration)
    resumed = plus_workflow.run_plus_flair_registration_checkpoint(
        FIXTURES / "fixed_image.nii.gz",
        FIXTURES / "moving_flair.nii.gz",
        output_dir,
        config=config,
        resume=True,
    )
    assert resumed.resumed
    np.testing.assert_array_equal(resumed.warped_flair.dataobj, first.warped_flair.dataobj)

    with pytest.raises(ValueError, match="inputs or settings differ"):
        plus_workflow.run_plus_flair_registration_checkpoint(
            FIXTURES / "fixed_image.nii.gz",
            FIXTURES / "moving_flair.nii.gz",
            output_dir,
            config=RegistrationConfig(transform_type="Rigid", metric="meansquares"),
            resume=True,
        )
