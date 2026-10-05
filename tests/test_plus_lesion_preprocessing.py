from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import save_nifti_ras
import pymacruise.plus_lesion_preprocessing as preprocessing
from pymacruise.n4 import N4Config
from pymacruise.registration import RegistrationConfig


FIXTURE = Path(__file__).parent / "data/synthetic/registration"


def test_plus_lesion_preprocessing_registers_corrects_strips_and_resumes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixed = nib.load(str(FIXTURE / "fixed_image.nii.gz"))
    assert isinstance(fixed, nib.Nifti1Image)
    brain = np.asarray(nib.load(str(FIXTURE / "fixed_mask.nii.gz")).dataobj) > 0
    mask_path = save_nifti_ras(
        nib.Nifti1Image(brain.astype(np.uint8), fixed.affine), tmp_path / "brain.nii.gz"
    )
    registration_config = RegistrationConfig(
        transform_type="Rigid",
        metric="mattes",
        sampling_rate=1.0,
        iterations=(20,),
        shrink_factors=(1,),
        smoothing_sigmas=(0,),
    )
    n4_config = N4Config(shrink_factor=1, iterations=(2,), spline_mesh=(1, 1, 1))
    output = tmp_path / "plus-lesion-prep"
    first = preprocessing.run_plus_lesion_preprocessing_checkpoint(
        FIXTURE / "fixed_image.nii.gz",
        FIXTURE / "moving_flair.nii.gz",
        mask_path,
        output,
        registration_config=registration_config,
        n4_config=n4_config,
    )
    assert not first.resumed
    assert all(path.is_file() for path in (first.t1_brain_path, first.flair_brain_path, first.modalities_path))
    t1 = nib.load(str(first.t1_brain_path))
    flair = nib.load(str(first.flair_brain_path))
    modalities = nib.load(str(first.modalities_path))
    assert t1.shape == flair.shape == fixed.shape
    assert modalities.shape == (*fixed.shape, 2)
    for image in (t1, flair, modalities):
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(image.dataobj)).all()
    np.testing.assert_array_equal(np.asarray(modalities.dataobj)[..., 0], np.asarray(t1.dataobj))
    np.testing.assert_array_equal(np.asarray(modalities.dataobj)[..., 1], np.asarray(flair.dataobj))
    np.testing.assert_array_equal(np.asarray(t1.dataobj)[~brain], 0.0)
    np.testing.assert_array_equal(np.asarray(flair.dataobj)[~brain], 0.0)
    manifest = json.loads(first.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["channel_order"] == ["T1", "FLAIR"]
    assert manifest["source_order"] == [
        "rigid FLAIR-to-T1 registration",
        "N4 T1",
        "N4 registered FLAIR",
        "brain-mask strip and negative cleanup",
    ]
    assert Path(manifest["registration_manifest"]["path"]).is_file()
    assert "NiftyReg/Aladin" in manifest["equivalence_note"]

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume should not repeat preprocessing")

    monkeypatch.setattr(preprocessing, "n4_bias_correct", unexpected)
    monkeypatch.setattr(preprocessing, "run_plus_flair_registration_checkpoint", unexpected)
    resumed = preprocessing.run_plus_lesion_preprocessing_checkpoint(
        FIXTURE / "fixed_image.nii.gz",
        FIXTURE / "moving_flair.nii.gz",
        mask_path,
        output,
        registration_config=registration_config,
        n4_config=n4_config,
        resume=True,
    )
    assert resumed.resumed
    assert resumed.modalities_path == first.modalities_path


def test_plus_lesion_preprocessing_rejects_mask_outside_t1_grid(tmp_path: Path) -> None:
    fixed = nib.load(str(FIXTURE / "fixed_image.nii.gz"))
    assert isinstance(fixed, nib.Nifti1Image)
    bad_mask = save_nifti_ras(
        nib.Nifti1Image(np.ones((16, 17, 17), dtype=np.uint8), fixed.affine),
        tmp_path / "bad-mask.nii.gz",
    )
    with pytest.raises(ValueError, match="brain mask must match the T1 reference grid"):
        preprocessing.run_plus_lesion_preprocessing_checkpoint(
            FIXTURE / "fixed_image.nii.gz",
            FIXTURE / "moving_flair.nii.gz",
            bad_mask,
            tmp_path / "unused",
            n4_config=N4Config(shrink_factor=1, iterations=(1,)),
        )


def test_plus_prepare_lesion_cli_routes_source_settings(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_prepare(t1, flair, mask, output, **kwargs):
        captured.update(t1=t1, flair=flair, mask=mask, output=output, **kwargs)
        return preprocessing.PlusLesionPreprocessingResult(
            tmp_path / "manifest.json",
            tmp_path / "t1.nii.gz",
            tmp_path / "flair.nii.gz",
            tmp_path / "modalities.nii.gz",
            resumed=False,
        )

    monkeypatch.setattr(preprocessing, "run_plus_lesion_preprocessing_checkpoint", fake_prepare)
    assert main(
        [
            "plus-prepare-lesion",
            "--t1", "t1.nii.gz",
            "--flair", "flair.nii.gz",
            "--brain-mask", "mask.nii.gz",
            "--output-dir", "out",
            "--registration-iterations", "5",
            "--shrink-factors", "1",
            "--smoothing-sigmas", "0",
            "--n4-shrink-factor", "1",
            "--n4-iterations", "2",
            "--n4-spline-mesh", "1", "1", "1",
            "--resume",
        ]
    ) == 0
    assert captured["registration_config"].iterations == (5,)
    assert captured["registration_config"].transform_type == "Rigid"
    assert captured["n4_config"].iterations == (2,)
    assert captured["n4_config"].spline_mesh == (1, 1, 1)
    assert captured["resume"] is True
