from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np


def _save(path: Path, values: np.ndarray, affine: np.ndarray) -> None:
    nib.save(nib.Nifti1Image(values, affine), path)


def test_module_cli_writes_class_member_products_as_ras_nifti(tmp_path: Path) -> None:
    shape = (5, 4, 3)
    # LAS fixture: the CLI warns/reorients the reference before writing RAS.
    affine = np.diag([-1.5, 2.0, 2.5, 1.0])
    affine[0, 3] = 6.0
    reference = tmp_path / "subject-t1.nii.gz"
    classification = tmp_path / "class.nii.gz"
    subtype = tmp_path / "subtype.nii.gz"
    membership = tmp_path / "membership.nii.gz"
    output_dir = tmp_path / "output"
    _save(reference, np.zeros(shape, dtype=np.float32), affine)
    _save(classification, np.arange(np.prod(shape), dtype=np.int8).reshape(shape), np.eye(4))
    _save(subtype, np.full(shape, 2, dtype=np.int16), np.eye(4))
    packed = np.arange(np.prod((*shape, 4)), dtype=np.float32).reshape((*shape, 4))
    _save(membership, packed, np.eye(4))

    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path(__file__).parents[1] / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.lesion_output_cli",
            "--reference",
            str(reference),
            "--classification",
            str(classification),
            "--lesion-subtype",
            str(subtype),
            "--mode",
            "class+member",
            "--mode-data",
            str(membership),
            "--input-name",
            "subject-t1.nii.gz",
            "--output-dir",
            str(output_dir),
        ],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "reorienting the image to RAS" in result.stderr
    saved = sorted(output_dir.glob("*.nii.gz"))
    assert {path.name for path in saved} == {
        "subject-t1_lToads_class.nii.gz",
        "subject-t1_lToads_lesions_hard.nii.gz",
        "subject-t1_lToads_mems.nii.gz",
    }
    written = nib.load(output_dir / "subject-t1_lToads_mems.nii.gz")
    assert written.shape == (*shape, 4)
    assert written.get_data_dtype() == np.dtype(np.float32)
    assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(written.dataobj, packed)


def test_module_cli_requires_mode_data_except_for_classification_only(tmp_path: Path) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path(__file__).parents[1] / "src")
    result = subprocess.run(
        [sys.executable, "-m", "pymacruise.lesion_output_cli", "--help"],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "already-computed LesionTOADS terminal arrays" in result.stdout
