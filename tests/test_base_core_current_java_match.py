"""Regression checks for the current saved base-run vs Java TGDM evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent / "data/synthetic/base_core_source_match"
MANIFEST = json.loads((ROOT / "comparison.json").read_text())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_source_match_fixture_hashes_and_ras_geometry() -> None:
    assert MANIFEST["synthetic_only"] is True
    assert MANIFEST["orientation"] == "RAS"
    for relative, expected_hash in MANIFEST["files"].items():
        path = ROOT / relative
        assert _sha256(path) == expected_hash

    for name in ("inner", "central", "outer"):
        java_img = nib.load(str(ROOT / "java_reference" / f"{name}.nii.gz"))
        python_img = nib.load(str(ROOT / "python_run" / f"{name}.nii.gz"))
        assert nib.aff2axcodes(java_img.affine) == ("R", "A", "S")
        assert nib.aff2axcodes(python_img.affine) == ("R", "A", "S")
        assert java_img.shape == python_img.shape == (172, 216, 164)
        np.testing.assert_allclose(java_img.affine, python_img.affine, rtol=0, atol=1e-5)


def test_current_run_tgdm_and_refined_label_metrics_meet_agreed_gates() -> None:
    for name in ("inner", "central", "outer"):
        metrics = MANIFEST["metrics"][name]
        assert metrics["dice"] >= 0.97
        assert metrics["mean_absolute_distance_mm"] <= 0.01
        assert metrics["hausdorff_mm"] <= 1.0

    refined = MANIFEST["metrics"]["refined_labels"]
    assert MANIFEST["refinement_source_faithful"] is True
    assert refined["changed_voxels"] == 7
    for label in ("1", "4", "100"):
        assert refined["labels"][label]["dice"] >= 0.97

