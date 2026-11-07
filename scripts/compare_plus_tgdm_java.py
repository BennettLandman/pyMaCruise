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

from pymacruise.surfaces import extract_levelset_surface
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
    peak_matches = re.findall(r"peak_sampled_rss_gb=([0-9.]+)", watchdog_log)
    if not peak_matches:
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
            java_surface = extract_levelset_surface(nib.Nifti1Image(java, affine))
            py_surface = extract_levelset_surface(nib.Nifti1Image(python, affine))
            java_vertices, java_faces, _, _ = marching_cubes(java, level=0.0)
            py_vertices, py_faces, _, _ = marching_cubes(python, level=0.0)
            item["dice_at_zero"] = dice_coefficient(java <= 0.0, python <= 0.0)
            item.update(mesh_surface_distance_metrics(
                java_surface.vertices, java_surface.faces,
                py_surface.vertices, py_surface.faces,
            ))
            item["marching_cubes_diagnostic"] = mesh_surface_distance_metrics(
                apply_affine(affine, java_vertices), java_faces,
                apply_affine(affine, py_vertices), py_faces,
            )
            item["thresholds_pass"] = (
                item["dice_at_zero"] >= 0.97
                and item["mean_absolute_distance_mm"] <= 0.01
                and item["hausdorff_mm"] <= 1.0
            )
            all_pass &= bool(item["thresholds_pass"])
        outputs[name] = item
    replay = provenance.get("python_tgdm_replay")
    input_mode = provenance.get("tgdm_input_mode", "initial_membership")
    if replay:
        comparison = "selected MaGenericTGDM Java null-bound wrapper sequence vs current Python null-bound TGDM replay on exact held-fixed inputs"
        scope = "TGDM only; saved upstream membership, tissue, ACE, skeleton, and GVF arrays are held fixed; no optional brain-mask level set is passed to either implementation"
    elif input_mode == "initial_membership":
        comparison = "selected MaGenericTGDM Java null-bound wrapper sequence vs Python null-bound replay from the same held-fixed initial membership"
        scope = "TGDM only; the same initial membership, tissue, ACE, skeleton, and GVF arrays are supplied to both implementations; each computes its selected initial signed-distance field; no optional brain-mask level set is passed"
    else:
        comparison = "selected MaGenericTGDM Java null-bound wrapper sequence vs saved Python TGDM outputs from a held-fixed signed initial level set"
        scope = "TGDM evolution only; the saved signed initial level set, tissue, ACE, skeleton, and GVF arrays are held fixed; no optional brain-mask level set is passed"
    report = {
        "synthetic_only": provenance["synthetic_only"],
        "evidence_category": provenance["evidence_category"],
        "comparison": comparison,
        "scope": scope,
        "tgdm_input_mode": input_mode,
        "python_tgdm_replay": replay,
        "tgdm_optional_brain_mask_level_set": provenance.get("tgdm_optional_brain_mask_level_set"),
        "brain_mask_role": provenance.get("brain_mask_role"),
        "replayed_saved_outputs_exactly": provenance.get("replayed_saved_outputs_exactly"),
        "saved_stage_difference_voxels": provenance.get("saved_stage_difference_voxels"),
        "held_fixed_input_binary_sha256": provenance.get("held_fixed_input_binary_sha256"),
        "workflow_run_signature": provenance["workflow_run_signature"],
        "full_workflow_run_signature": provenance["full_workflow_run_signature"],
        "implementation_signature": provenance["implementation_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "surface_extractor": "source-matched JIST IsoSurfaceOnGrid connectivity-18 port",
        "surface_gate_note": "Primary distance metrics use production JIST-compatible surfaces; marching_cubes_diagnostic is informational only.",
        "java_peak_sampled_rss_gb": float(peak_matches[-1]),
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
