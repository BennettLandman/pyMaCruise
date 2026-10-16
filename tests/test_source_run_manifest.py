from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_source_reference_run.py"


def _run_audit(run_dir: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(run_dir), *args],
        check=False,
        capture_output=True,
        text=True,
    )


def test_source_run_manifest_records_hashes_and_ras_geometry(tmp_path: Path):
    run_dir = tmp_path / "source-run"
    run_dir.mkdir()
    affine = np.eye(4)
    values = np.arange(24, dtype=np.float32).reshape(3, 4, 2)
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, str(run_dir / "inner.nii.gz"))
    (run_dir / "RUN_INFO.md").write_text("source revision: example\n", encoding="utf-8")
    (run_dir / "inner.vtk").write_text("synthetic surface placeholder\n", encoding="utf-8")

    result = _run_audit(run_dir, "--require-ras", "--require-run-info")
    assert result.returncode == 0, result.stderr
    manifest = json.loads((run_dir / "source-run-manifest.json").read_text())
    assert manifest["passed"] is True
    assert manifest["file_count"] == 3
    assert manifest["nifti_count"] == 1
    assert manifest["has_run_info"] is True
    image_record = manifest["files"]["inner.nii.gz"]["nifti"]
    assert image_record["shape"] == [3, 4, 2]
    assert image_record["axis_codes"] == ["R", "A", "S"]
    assert image_record["finite_voxels"] is True
    assert image_record["qform_code"] == 1
    assert image_record["sform_code"] == 1


def test_source_run_manifest_rejects_nonras_and_uncompressed_images(tmp_path: Path):
    run_dir = tmp_path / "source-run"
    run_dir.mkdir()
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    values = np.ones((3, 4, 2), dtype=np.float32)
    image = nib.Nifti1Image(values, affine)
    nib.save(image, str(run_dir / "left-right.nii.gz"))
    nib.save(image, str(run_dir / "uncompressed.nii"))

    result = _run_audit(run_dir, "--require-ras")
    assert result.returncode == 1
    manifest = json.loads((run_dir / "source-run-manifest.json").read_text())
    assert manifest["passed"] is False
    assert any("non-RAS NIfTI" in error for error in manifest["errors"])
    assert any("uncompressed NIfTI" in error for error in manifest["errors"])


def test_source_run_manifest_reports_nonfinite_voxels_without_nonstandard_json(tmp_path: Path):
    run_dir = tmp_path / "source-run"
    run_dir.mkdir()
    values = np.array([[[1.0], [np.nan]], [[np.inf], [-np.inf]]], dtype=np.float32)
    nib.save(nib.Nifti1Image(values, np.eye(4)), str(run_dir / "nonfinite.nii.gz"))

    result = _run_audit(run_dir)
    assert result.returncode == 1
    manifest = json.loads((run_dir / "source-run-manifest.json").read_text())
    image_record = manifest["files"]["nonfinite.nii.gz"]["nifti"]
    assert image_record["finite_voxels"] is False
    assert image_record["minimum"] == 1.0
    assert image_record["maximum"] == 1.0
    assert "nonfinite NIfTI voxel data" in manifest["errors"][0]


def test_source_run_manifest_rejects_complex_nifti(tmp_path: Path):
    run_dir = tmp_path / "source-run"
    run_dir.mkdir()
    values = np.ones((2, 2, 2), dtype=np.complex64) * (1 + 2j)
    nib.save(nib.Nifti1Image(values, np.eye(4)), str(run_dir / "complex.nii.gz"))

    result = _run_audit(run_dir)
    assert result.returncode == 1
    manifest = json.loads((run_dir / "source-run-manifest.json").read_text())
    assert "could not audit NIfTI complex.nii.gz" in manifest["errors"][0]
    assert "nonnumeric NIfTI datatype" in manifest["errors"][0]
