#!/usr/bin/env python3
"""Compare the synthetic base MaCRUISE TGDM core with selected Java source calls."""

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

from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_images  # noqa: E402
from pymacruise.validation import (  # noqa: E402
    dice_coefficient,
    mesh_surface_distance_metrics,
    voxel_error_metrics,
)

DEFAULT_INPUT = REPO_ROOT / "tests/data/synthetic/tgdm_core_brain_resolution/java_inputs"
DEFAULT_REFERENCE = REPO_ROOT / "tests/data/synthetic/tgdm_core_brain_resolution/java_reference"


def world_mesh(values: np.ndarray, affine: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vertices, faces, _, _ = marching_cubes(values, level=0.0)
    vertices = nib.affines.apply_affine(affine, vertices)
    return vertices, faces


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--reference-root", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--candidate-root", type=Path,
                        help="cache Python stage arrays here, or load existing cached arrays")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    inputs = {name: nib.load(str(args.input_root / f"{name}.nii.gz")) for name in (
        "initial_membership", "white_matter", "gray_matter", "ace_gray_matter", "ace_skeleton", "gvf"
    )}
    grid = inputs["initial_membership"]
    params = TGDMBaselineParameters(
        white_matter_iso=0.6,
        initial_iterations=1,
        inner_iterations=1,
        central_shrink_iterations=1,
        central_iterations=1,
        outer_iterations=1,
    )
    result = tgdm_nested_images(
        inputs["initial_membership"], inputs["white_matter"], inputs["gray_matter"],
        inputs["ace_gray_matter"], ace_skeleton=inputs["ace_skeleton"],
        gvf=inputs["gvf"], parameters=params,
    )
    stages = {
        "signed_seed": result.execution.initial_level_set,
        "smoothed_initial": result.execution.smoothed_initial.level_set,
        "initial_inner": result.execution.initial_inner.level_set,
        "tissue_sign": result.tissue_sign,
        "central": result.execution.central.level_set,
        "central_shrink": result.execution.central_shrink.level_set,
        "inner": result.execution.inner.level_set,
        "outer": result.execution.outer.level_set,
    }
    if args.candidate_root:
        args.candidate_root.mkdir(parents=True, exist_ok=True)
        for name, values in stages.items():
            path = args.candidate_root / f"{name}.nii.gz"
            try:
                cached = np.asarray(nib.load(str(path)).dataobj) if path.exists() else None
            except (OSError, ValueError):
                cached = None
            if cached is not None:
                stages[name] = cached
            else:
                if hasattr(values, "dataobj"):
                    values = np.asarray(values.dataobj)
                else:
                    values = np.asarray(values)
                output = nib.Nifti1Image(values, grid.affine, grid.header.copy())
                nib.save(output, str(path))
    report: dict[str, object] = {
        "synthetic_only": True,
        "comparison": "Python selected TGDM core vs supplied MaGenericTGDM/DistanceField Java calls",
        "orientation": "RAS",
        "shape": list(grid.shape),
        "voxel_spacing_mm": [float(x) for x in nib.affines.voxel_sizes(grid.affine)],
        "parameters": {
            "white_matter_iso": 0.6,
            "gray_matter_iso": 0.4,
            "initial_iterations": 1,
            "inner_iterations": 1,
            "central_iterations": 1,
            "outer_wrapper_iterations": 1,
            "nonzero_gvf": True,
        },
        "acceptance": {"mask_dice_minimum": 0.97, "symmetric_mean_surface_distance_mm_maximum": 0.01,
                       "symmetric_hausdorff_mm_maximum": 1.0},
        "stages": {},
    }
    all_pass = True
    for name, candidate in stages.items():
        reference_image = nib.load(str(args.reference_root / f"{name}.nii.gz"))
        if nib.aff2axcodes(reference_image.affine) != ("R", "A", "S"):
            raise ValueError(f"Java reference is not RAS: {name}")
        reference = np.asarray(reference_image.dataobj)
        if name == "tissue_sign":
            candidate_array = np.asarray(candidate.dataobj if hasattr(candidate, "dataobj") else candidate)
            metrics = {"exact": bool(np.array_equal(candidate_array, reference)),
                       **voxel_error_metrics(candidate_array, reference)}
            all_pass &= metrics["exact"]
        else:
            candidate = np.asarray(candidate, dtype=np.float32)
            reference = np.asarray(reference, dtype=np.float32)
            metrics = {
                "dice": dice_coefficient(candidate <= 0.0, reference <= 0.0),
                **voxel_error_metrics(candidate, reference),
            }
            # Full triangle-surface distances are reserved for the three
            # delivered cortical surfaces; intermediate fields still receive
            # full-grid voxel errors and sign-mask Dice without repeating a
            # costly mesh comparison at each evolution step.
            if name in {"central", "inner", "outer"}:
                candidate_vertices, candidate_faces = world_mesh(candidate, grid.affine)
                reference_vertices, reference_faces = world_mesh(reference, grid.affine)
                metrics.update(mesh_surface_distance_metrics(
                    candidate_vertices, candidate_faces,
                    reference_vertices, reference_faces,
                ))
                all_pass &= (
                    metrics["dice"] >= 0.97
                    and metrics["mean_absolute_distance_mm"] <= 0.01
                    and metrics["hausdorff_mm"] <= 1.0
                )
            else:
                all_pass &= metrics["dice"] >= 0.97
        report["stages"][name] = metrics
    report["passed"] = bool(all_pass)
    serialized = json.dumps(report, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n")
    if not all_pass:
        raise SystemExit("TGDM core Java source comparison did not meet acceptance thresholds")


if __name__ == "__main__":
    main()
