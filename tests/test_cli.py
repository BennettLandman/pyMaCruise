from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import main


def test_cli_version_and_prepare_raw_commands(tmp_path, capsys) -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert "pymacruise" in capsys.readouterr().out

    affine = np.eye(4)
    t1_path = tmp_path / "t1.nii.gz"
    maseg_path = tmp_path / "maseg.nii.gz"
    mask_path = tmp_path / "brain-mask.nii.gz"
    nib.save(nib.Nifti1Image(np.ones((4, 5, 6), dtype=np.float32), affine), t1_path)
    labels = np.zeros((4, 5, 6), dtype=np.uint8)
    labels[1:3, 1:4, 2:5] = 100
    nib.save(nib.Nifti1Image(labels, affine), maseg_path)
    brain_mask = np.zeros((4, 5, 6), dtype=np.uint8)
    brain_mask[1:3, 1:4, 1:5] = 255
    nib.save(nib.Nifti1Image(brain_mask, affine), mask_path)

    assert main(
        [
            "prepare-raw",
            "--t1",
            str(t1_path),
            "--maseg",
            str(maseg_path),
            "--brain-mask",
            str(mask_path),
            "--output-dir",
            str(tmp_path / "prepared"),
        ]
    ) == 0
    assert (tmp_path / "prepared" / "T1.nii.gz").is_file()
    assert (tmp_path / "prepared" / "Seg.nii.gz").is_file()
    assert (tmp_path / "prepared" / "BrainMask.nii.gz").is_file()
