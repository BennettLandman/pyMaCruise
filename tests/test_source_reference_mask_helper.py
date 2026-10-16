from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "source-reference-run/matlab/apply_brain_mask_nifti.py"


def test_source_reference_mask_helper_preserves_values_and_ras_geometry(tmp_path: Path):
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    affine[:3, 3] = (-1.0, -2.0, -3.0)
    values = np.arange(24, dtype=np.float32).reshape(3, 4, 2)
    mask_values = np.zeros(values.shape, dtype=np.uint8)
    mask_values[1:, 1:3, :] = 1
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    mask = nib.Nifti1Image(mask_values, affine)
    mask.set_qform(affine, code=1)
    mask.set_sform(affine, code=1)
    image_path, mask_path = tmp_path / "n4.nii.gz", tmp_path / "mask.nii.gz"
    output_path = tmp_path / "masked.nii.gz"
    nib.save(image, str(image_path))
    nib.save(mask, str(mask_path))

    subprocess.run(
        [
            sys.executable,
            str(HELPER),
            "--image",
            str(image_path),
            "--mask",
            str(mask_path),
            "--output",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    output = nib.load(str(output_path))
    actual = np.asarray(output.dataobj)
    assert np.array_equal(actual[mask_values == 0], np.zeros(np.count_nonzero(mask_values == 0)))
    assert np.array_equal(actual[mask_values != 0], values[mask_values != 0])
    assert output.shape == values.shape
    assert np.allclose(output.affine, affine, rtol=0, atol=1e-5)
    assert nib.aff2axcodes(output.affine) == ("R", "A", "S")
    _, qcode = output.get_qform(coded=True)
    _, scode = output.get_sform(coded=True)
    assert qcode == 1 and scode == 1


def test_source_reference_mask_helper_rejects_mismatched_grid(tmp_path: Path):
    affine = np.eye(4)
    image_path, mask_path = tmp_path / "n4.nii.gz", tmp_path / "mask.nii.gz"
    nib.save(nib.Nifti1Image(np.ones((3, 4, 2), dtype=np.float32), affine), str(image_path))
    shifted = affine.copy()
    shifted[0, 3] = 2.0
    nib.save(nib.Nifti1Image(np.ones((3, 4, 2), dtype=np.uint8), shifted), str(mask_path))
    result = subprocess.run(
        [
            sys.executable,
            str(HELPER),
            "--image",
            str(image_path),
            "--mask",
            str(mask_path),
            "--output",
            str(tmp_path / "masked.nii.gz"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "mask shape and affine must match" in result.stderr
