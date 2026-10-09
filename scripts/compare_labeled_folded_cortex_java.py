#!/usr/bin/env python3
"""Compare full-profile Python TGDM with the selected Java harness outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics, voxel_error_metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_root", type=Path)
    parser.add_argument("java_root", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    stages = (
        "signed_seed", "smoothed_initial", "initial_inner", "tissue_sign",
        "central", "central_shrink", "inner", "outer",
    )
    report: dict[str, object] = {
        "evidence_category": "selected Java/Python source comparison on one synthetic labeled folded-cortex case",
        "source_profile": "full TGDMBaselineParameters with run_ToadsCruise WM iso 0.6 and explicit signed brain bound",
        "orientation": "RAS",
        "stages": {},
        "acceptance": {
            "final_surface_dice_minimum": 0.97,
            "symmetric_mean_surface_distance_mm_maximum": 0.01,
            "symmetric_hausdorff_mm_maximum": 1.0,
        },
    }
    stage_reports = {}
    final_pass = True
    for name in stages:
        py_image = nib.load(str(args.python_root / f"{name}.nii.gz"))
        java_image = nib.load(str(args.java_root / f"{name}.nii.gz"))
        if nib.aff2axcodes(py_image.affine) != ("R", "A", "S") or nib.aff2axcodes(java_image.affine) != ("R", "A", "S"):
            raise ValueError(f"{name}: both images must be RAS")
        if py_image.shape != java_image.shape or not np.allclose(py_image.affine, java_image.affine):
            raise ValueError(f"{name}: Python and Java grids differ")
        py = np.asarray(py_image.dataobj, dtype=np.float32)
        java = np.asarray(java_image.dataobj, dtype=np.float32)
        record: dict[str, object] = {
            "python_sha256": __import__("hashlib").sha256((args.python_root / f"{name}.nii.gz").read_bytes()).hexdigest(),
            "java_sha256": __import__("hashlib").sha256((args.java_root / f"{name}.nii.gz").read_bytes()).hexdigest(),
            **voxel_error_metrics(py, java),
        }
        if name != "tissue_sign":
            record["region_dice"] = dice_coefficient(py <= 0, java <= 0)
        if name in {"initial_inner", "central", "inner", "outer"}:
            pv, pf, _, _ = marching_cubes(py, level=0.0)
            jv, jf, _, _ = marching_cubes(java, level=0.0)
            pv = nib.affines.apply_affine(py_image.affine, pv)
            jv = nib.affines.apply_affine(java_image.affine, jv)
            record["surface"] = mesh_surface_distance_metrics(pv, pf, jv, jf)
        if name in {"inner", "central", "outer"}:
            surface = record["surface"]
            passed = (
                record["region_dice"] >= report["acceptance"]["final_surface_dice_minimum"]
                and surface["mean_absolute_distance_mm"]
                <= report["acceptance"]["symmetric_mean_surface_distance_mm_maximum"]
                and surface["hausdorff_mm"] <= report["acceptance"]["symmetric_hausdorff_mm_maximum"]
            )
            record["accepted"] = bool(passed)
            final_pass &= bool(passed)
        stage_reports[name] = record
    report["stages"] = stage_reports
    report["final_surfaces_passed"] = final_pass
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if not final_pass:
        raise SystemExit("one or more final Java/Python surface gates failed")


if __name__ == "__main__":
    main()
