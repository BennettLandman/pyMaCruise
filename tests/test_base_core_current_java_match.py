"""Regression checks for the current saved base-run vs Java TGDM evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

from pymacruise.validation import (
    dice_coefficient,
    mesh_surface_distance_metrics,
)


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
        java_img = nib.load(str(ROOT / "java_reference" / f"{name}.nii.gz"))
        python_img = nib.load(str(ROOT / "python_run" / f"{name}.nii.gz"))
        java_values = np.asarray(java_img.dataobj, dtype=np.float32)
        python_values = np.asarray(python_img.dataobj, dtype=np.float32)
        java_vertices, java_faces, _, _ = marching_cubes(java_values, level=0.0)
        python_vertices, python_faces, _, _ = marching_cubes(python_values, level=0.0)
        java_vertices = nib.affines.apply_affine(java_img.affine, java_vertices)
        python_vertices = nib.affines.apply_affine(python_img.affine, python_vertices)
        measured = {
            "dice": dice_coefficient(python_values <= 0, java_values <= 0),
            **mesh_surface_distance_metrics(
                python_vertices, python_faces, java_vertices, java_faces
            ),
        }
        recorded = MANIFEST["metrics"][name]
        for key in ("dice", "mean_absolute_distance_mm", "hausdorff_mm"):
            np.testing.assert_allclose(measured[key], recorded[key], rtol=1e-5, atol=1e-7)
        assert measured["dice"] >= 0.97
        assert measured["mean_absolute_distance_mm"] <= 0.01
        assert measured["hausdorff_mm"] <= 1.0

    refined = MANIFEST["metrics"]["refined_labels"]
    assert MANIFEST["refinement_source_faithful"] is True
    assert refined["changed_voxels"] == 7
    for label in ("1", "4", "100"):
        assert refined["labels"][label]["dice"] >= 0.97
