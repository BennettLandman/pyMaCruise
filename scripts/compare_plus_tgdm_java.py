#!/usr/bin/env python3
"""Compare held-fixed plus TGDM outputs with the selected Java wrapper sequence."""
from __future__ import annotations
import argparse
import gzip
import json
import re
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.affines import apply_affine
from skimage.measure import marching_cubes

from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics

NAMES = ("initial_level_set", "smoothed_initial", "initial_inner", "central", "inner", "outer")
SURFACES = ("smoothed_initial", "initial_inner", "central", "inner", "outer")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    workflow = args.workflow_root.resolve()
    stage = workflow / "stages/topology_ace_gvf_tgdm" / json.loads((workflow / "run.json").read_text())["stages"]["topology_ace_gvf_tgdm"]["attempt"]
    provenance = json.loads((args.reference_dir / "input-provenance.json").read_text())
    shape = tuple(provenance["shape_xyz"])
    affine = np.asarray(provenance["affine_ras"], dtype=np.float64)
    watchdog_log = (args.reference_dir / "java-watchdog.log").read_text()
    peak_match = re.search(r"peak_sampled_rss_gb=([0-9.]+)", watchdog_log)
    if peak_match is None:
        raise SystemExit("Java watchdog log does not contain a peak RSS measurement")
    outputs = {}
    all_pass = True
    for name in NAMES:
        with gzip.open(args.reference_dir / f"java/{name}.raw.gz", "rb") as stream:
            java = np.frombuffer(stream.read(), dtype=">f4").astype(np.float32).reshape(shape)
        python_path = args.reference_dir / (
            "python-smoothed_initial.npy" if name == "smoothed_initial" else f"python-{name}.npy"
        )
        python = np.load(python_path).astype(np.float32)
        diff = np.abs(java.astype(np.float64) - python.astype(np.float64))
        item = {
            "exact": bool(np.array_equal(java, python)),
            "changed_voxels": int(np.count_nonzero(java != python)),
            "voxel_count": int(java.size),
            "mean_absolute_voxel_error": float(diff.mean()),
            "maximum_absolute_voxel_error": float(diff.max(initial=0.0)),
            "sign_changed_voxels": int(np.count_nonzero(np.signbit(java) != np.signbit(python))),
        }
        if name in SURFACES:
            java_vertices, java_faces, _, _ = marching_cubes(java, level=0.0)
            py_vertices, py_faces, _, _ = marching_cubes(python, level=0.0)
            item["dice_at_zero"] = dice_coefficient(java <= 0.0, python <= 0.0)
            item.update(mesh_surface_distance_metrics(
                apply_affine(affine, java_vertices), java_faces,
                apply_affine(affine, py_vertices), py_faces,
            ))
            item["thresholds_pass"] = (
                item["dice_at_zero"] >= 0.97
                and item["mean_absolute_distance_mm"] <= 0.01
                and item["hausdorff_mm"] <= 1.0
            )
            all_pass &= bool(item["thresholds_pass"])
        outputs[name] = item
    report = {
        "synthetic_only": True,
        "comparison": "selected MaGenericTGDM Java plus-wrapper sequence vs saved Python plus workflow on exact held-fixed inputs",
        "scope": "TGDM only; the source signed-distance seed is reused from the hash-checked Python stage, tissue/ACE/GVF inputs are reconstructed on the same crop grid, and TOADS/LesionTOADS remain held fixed",
        "tgdm_input_mode": provenance.get("tgdm_input_mode", "initial_membership"),
        "workflow_run_signature": provenance["workflow_run_signature"],
        "full_workflow_run_signature": provenance["full_workflow_run_signature"],
        "implementation_signature": provenance["implementation_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "java_peak_sampled_rss_gb": float(peak_match.group(1)),
        "java_rss_limit_gb": 8,
        "affine_ras": affine.tolist(),
        "acceptance": {"dice_minimum": 0.97, "symmetric_mean_surface_distance_mm_maximum": 0.01, "symmetric_hausdorff_mm_maximum": 1.0},
        "outputs": outputs,
        "surface_thresholds_pass": bool(all_pass),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    if not all_pass:
        raise SystemExit("one or more held-fixed TGDM surfaces missed an agreed tolerance")


if __name__ == "__main__":
    main()
