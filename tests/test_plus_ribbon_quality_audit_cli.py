from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _record(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _write_labeled_sphere_case(folder: Path) -> Path:
    shape = (41, 41, 41)
    grid = np.indices(shape, dtype=np.float32)
    radius = np.sqrt(sum((axis - np.float32(20.0)) ** 2 for axis in grid))
    affine = np.eye(4, dtype=np.float32)
    levels = {
        "plus_inner_level_set": radius - np.float32(9.0),
        "plus_central_level_set": radius - np.float32(11.5),
        "plus_outer_level_set": radius - np.float32(14.0),
    }
    outputs = {}
    for name, values in levels.items():
        path = folder / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(values.astype(np.float32), affine), path)
        outputs[name] = _record(path)

    # The surfaces have a clean 5 mm nested shell, but intentionally put GM
    # outside the outer surface so the label/support criterion must fail.
    labels = np.zeros(shape, dtype=np.uint8)
    labels[radius < 8.0] = 44
    labels[(radius >= 15.0) & (radius < 16.0)] = 5
    labels[radius >= 16.0] = 100
    label_path = folder / "maseg.nii.gz"
    nib.save(nib.Nifti1Image(labels, affine), label_path)
    manifest = {
        "status": "complete",
        "stage": "plus-full-workflow",
        "terminal": True,
        "run_signature": "synthetic-failing-ribbon-quality-case",
        "inputs": {"maseg": _record(label_path)},
        "final_outputs": outputs,
    }
    manifest_path = folder / "plus-full-workflow.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return manifest_path


def test_require_quality_writes_failed_audit_and_returns_nonzero(tmp_path: Path) -> None:
    manifest = _write_labeled_sphere_case(tmp_path)
    report_path = tmp_path / "quality.json"

    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "audit_plus_ribbon_quality.py"),
            str(manifest),
            "--report",
            str(report_path),
            "--require-quality",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2, completed.stdout + completed.stderr
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["passed"] is False
    assert report["checks"]["surfaces_closed_connected_euler_two"] is True
    assert report["checks"]["nested_voxel_regions"] is True
    assert report["checks"]["labeled_tissue_regions_agree_with_surfaces"] is False
    assert "ribbon-quality requirements failed" in completed.stderr
