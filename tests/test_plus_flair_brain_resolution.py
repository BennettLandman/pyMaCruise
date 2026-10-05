from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.plus_workflow import run_plus_flair_registration_checkpoint
from pymacruise.registration import RegistrationConfig


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "plus_flair_brain_resolution"


def test_brain_resolution_synthetic_flair_checkpoint_registers_and_persists(tmp_path: Path) -> None:
    t1_path = FIXTURE / "t1.nii.gz"
    flair_path = FIXTURE / "flair_shifted.nii.gz"
    config = RegistrationConfig(
        metric="meansquares",
        sampling_rate=0.35,
        iterations=(40, 20, 0),
        shrink_factors=(4, 2, 1),
        smoothing_sigmas=(2.0, 1.0, 0.0),
        deterministic=True,
        single_precision=True,
    )
    result = run_plus_flair_registration_checkpoint(
        t1_path,
        flair_path,
        tmp_path / "flair-checkpoint",
        config=config,
    )
    fixed = np.asarray(nib.load(str(t1_path)).dataobj, dtype=np.float32)
    moving = np.asarray(nib.load(str(flair_path)).dataobj, dtype=np.float32)
    warped = np.asarray(result.warped_flair.dataobj, dtype=np.float32)
    before = float(np.corrcoef(fixed.ravel(), moving.ravel())[0, 1])
    after = float(np.corrcoef(fixed.ravel(), warped.ravel())[0, 1])

    assert result.warped_flair.shape == fixed.shape == (172, 216, 164)
    assert nib.aff2axcodes(result.warped_flair.affine) == ("R", "A", "S")
    assert np.isfinite(warped).all()
    assert after > before + 0.03
    assert result.warped_flair_path.is_file()
    assert result.manifest_path.is_file()
    assert result.forward_transforms and result.inverse_transforms
    assert all(path.is_file() for path in (*result.forward_transforms, *result.inverse_transforms))
