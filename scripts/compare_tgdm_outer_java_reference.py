#!/usr/bin/env python3
"""Compare the selected Python outer solver with Java from a shared central field."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from pymacruise.tgdm import (  # noqa: E402
    TGDMBaselineParameters,
    TGDMProfileSettings,
    tgdm_evolve,
    tgdm_outer_pressure_field,
)
from pymacruise.validation import (  # noqa: E402
    dice_coefficient,
    mesh_surface_distance_metrics,
    voxel_error_metrics,
)

DATA_ROOT = REPO_ROOT / "tests/data/synthetic/tgdm_core_brain_resolution"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DATA_ROOT / "java_inputs")
    parser.add_argument("--java-core-root", type=Path, default=DATA_ROOT / "java_reference")
    parser.add_argument("--java-outer-root", type=Path, default=DATA_ROOT / "java_outer_isolation")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    central_image = nib.load(str(args.java_core_root / "central.nii.gz"))
    gm_image = nib.load(str(args.input_root / "gray_matter.nii.gz"))
    skeleton_image = nib.load(str(args.input_root / "ace_skeleton.nii.gz"))
    java_output = nib.load(str(args.java_outer_root / "outer_from_java_central.nii.gz"))
    java_pre_reinit = nib.load(str(args.java_outer_root / "outer_before_reinit.nii.gz"))
    for label, image in (("central", central_image), ("GM", gm_image),
                         ("skeleton", skeleton_image), ("Java output", java_output),
                         ("Java pre-reinitialization output", java_pre_reinit)):
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{label} must be RAS")
        if image.shape != central_image.shape or not np.allclose(image.affine, central_image.affine):
            raise ValueError(f"{label} must share the central RAS grid")

    central = np.asarray(central_image.dataobj, dtype=np.float32)
    gray = np.asarray(gm_image.dataobj, dtype=np.float32)
    skeleton = np.asarray(skeleton_image.dataobj)
    settings = TGDMBaselineParameters()
    f32 = lambda value: float(np.float32(value))
    # OuterCortexProfile passes 15f, but ProfileTGDM's selected constructor
    # never assigns the finalMaxSignDist argument. Java therefore uses 20.
    outer_profile = TGDMProfileSettings(
        f32(4.7), f32(3.6), f32(np.float32(4.7) + np.float32(1.5)),
        f32(20.0), 30, 30, f32(0.1), f32(0.045),
    )
    force = tgdm_outer_pressure_field(
        central, gray, skeleton, gray_matter_iso=settings.gray_matter_iso
    )
    python_evolution = tgdm_evolve(
        central,
        profile="outer",
        curvature_force=-settings.outer_curvature,
        pressure_force=settings.outer_pressure,
        external_force=settings.outer_external,
        region_force=force,
        inner_bound=central,
        skeletal_bound=skeleton,
        outer_iterations=1,
        settings=outer_profile,
        record_passes=True,
    )
    python = python_evolution.level_set
    reference = np.asarray(java_output.dataobj, dtype=np.float32)
    python_pre = python_evolution.pre_reinitialization_level_sets[0]
    java_pre = np.asarray(java_pre_reinit.dataobj, dtype=np.float32)
    debug_root = REPO_ROOT / "build/tgdm-outer-isolation"
    debug_root.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(python, central_image.affine, central_image.header.copy()),
             str(debug_root / "outer_python.nii.gz"))
    nib.save(nib.Nifti1Image(python_pre, central_image.affine, central_image.header.copy()),
             str(debug_root / "outer_python_before_reinit.nii.gz"))
    candidate_vertices, candidate_faces, _, _ = marching_cubes(python, level=0.0)
    reference_vertices, reference_faces, _, _ = marching_cubes(reference, level=0.0)
    candidate_vertices = nib.affines.apply_affine(central_image.affine, candidate_vertices)
    reference_vertices = nib.affines.apply_affine(java_output.affine, reference_vertices)
    metrics = {
        "dice": dice_coefficient(python <= 0, reference <= 0),
        **voxel_error_metrics(python, reference),
        **mesh_surface_distance_metrics(candidate_vertices, candidate_faces,
                                        reference_vertices, reference_faces),
    }
    pre_reinit_metrics = voxel_error_metrics(python_pre, java_pre)
    pre_difference = np.abs(python_pre - java_pre)
    pre_outlier = tuple(int(value) for value in np.unravel_index(np.argmax(pre_difference), pre_difference.shape))
    pre_sign_disagreements = int(np.count_nonzero((python_pre <= 0) != (java_pre <= 0)))
    report = {
        "synthetic_only": True,
        "comparison": "isolated selected outer solver; identical Java central level set, GM, skeleton and source profile",
        "java_profile_note": "OuterCortexProfile passes 15f, but selected ProfileTGDM constructor does not assign finalMaxSignDist; the active field remains 20, confirmed by DistanceField runtime output",
        "orientation": "RAS",
        "shape": list(central_image.shape),
        "metrics": metrics,
        "pre_reinitialization_metrics": pre_reinit_metrics,
        "pre_reinitialization_max_difference_voxel": list(pre_outlier),
        "pre_reinitialization_python_value": float(python_pre[pre_outlier]),
        "pre_reinitialization_java_value": float(java_pre[pre_outlier]),
        "pre_reinitialization_sign_disagreements": pre_sign_disagreements,
        "acceptance": {"dice_minimum": 0.97, "symmetric_mean_surface_distance_mm_maximum": 0.01,
                       "symmetric_hausdorff_mm_maximum": 1.0},
        "passed": bool(metrics["dice"] >= 0.97
                       and metrics["mean_absolute_distance_mm"] <= 0.01
                       and metrics["hausdorff_mm"] <= 1.0),
    }
    serialized = json.dumps(report, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n")
    if not report["passed"]:
        raise SystemExit("Isolated Python outer-surface solver differs from Java reference")


if __name__ == "__main__":
    main()
