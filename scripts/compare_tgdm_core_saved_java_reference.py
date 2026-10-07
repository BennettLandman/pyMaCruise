#!/usr/bin/env python3
"""Compare saved synthetic base-run cortical outputs with Java TGDM references."""

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

from pymacruise.cropping import CropParameters, uncrop_image  # noqa: E402
from pymacruise.validation import (  # noqa: E402
    dice_coefficient,
    mesh_surface_distance_metrics,
    voxel_error_metrics,
)

DEFAULT_CORE = REPO_ROOT / "build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-02"
DEFAULT_REFERENCE = REPO_ROOT / "tests/data/synthetic/tgdm_core_brain_resolution/java_reference"
OUTPUTS = ("inner", "central", "outer")


def world_mesh(values: np.ndarray, affine: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vertices, faces, _, _ = marching_cubes(values, level=0.0)
    return nib.affines.apply_affine(affine, vertices), faces


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, default=DEFAULT_CORE,
                        help="saved output of the restartable synthetic base-run cruise-core stage")
    parser.add_argument("--reference-root", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    provenance = json.loads((args.core_dir / "provenance.json").read_text())
    source_prov = provenance.get("parameters", {})
    crop_bounds = tuple(tuple(int(v) for v in pair) for pair in source_prov["crop_bounds"])
    prototype = nib.load(str(args.core_dir / "inner.nii.gz"))
    crop_parameters = CropParameters(
        original_shape=tuple(int(value) for value in prototype.shape),
        bounds=crop_bounds,
        background=float(source_prov["uncrop_background"]),
        original_affine=prototype.affine.copy(),
    )
    if nib.aff2axcodes(prototype.affine) != ("R", "A", "S"):
        raise ValueError("saved Python base-run output must be RAS")

    stages: dict[str, object] = {}
    passed = True
    for name in OUTPUTS:
        python_image = nib.load(str(args.core_dir / f"{name}.nii.gz"))
        java_crop = nib.load(str(args.reference_root / f"{name}.nii.gz"))
        java_native = uncrop_image(java_crop, crop_parameters)
        if python_image.shape != java_native.shape or not np.allclose(
            python_image.affine, java_native.affine, rtol=0.0, atol=1e-5
        ):
            raise ValueError(f"{name}: Python and Java outputs do not share the native RAS grid")
        candidate = np.asarray(python_image.dataobj, dtype=np.float32)
        reference = np.asarray(java_native.dataobj, dtype=np.float32)
        candidate_vertices, candidate_faces = world_mesh(candidate, python_image.affine)
        reference_vertices, reference_faces = world_mesh(reference, java_native.affine)
        metrics = {
            "dice": dice_coefficient(candidate <= 0.0, reference <= 0.0),
            **voxel_error_metrics(candidate, reference),
            **mesh_surface_distance_metrics(candidate_vertices, candidate_faces,
                                            reference_vertices, reference_faces),
        }
        passed &= (
            metrics["dice"] >= 0.97
            and metrics["mean_absolute_distance_mm"] <= 0.01
            and metrics["hausdorff_mm"] <= 1.0
        )
        stages[name] = metrics

    report = {
        "synthetic_only": True,
        "comparison": "saved Python synthetic base-run core outputs vs selected MaGenericTGDM Java outputs",
        "source_harness": "TGDMCoreReferenceHarness; same RAS crop inputs and one-iteration wrapper profile",
        "python_source": str(args.core_dir),
        "orientation": "RAS",
        "shape": list(prototype.shape),
        "voxel_spacing_mm": [float(value) for value in nib.affines.voxel_sizes(prototype.affine)],
        "acceptance": {"mask_dice_minimum": 0.97,
                       "symmetric_mean_surface_distance_mm_maximum": 0.01,
                       "symmetric_hausdorff_mm_maximum": 1.0},
        "stages": stages,
        "passed": bool(passed),
        "scope_note": "source comparison covers delivered inner, central and outer masks on native RAS grid; intermediate evolution fields are retained in the Java fixture but not re-run here",
    }
    serialized = json.dumps(report, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n")
    if not passed:
        raise SystemExit("Saved synthetic TGDM core outputs do not meet the numerical acceptance gates")


if __name__ == "__main__":
    main()
