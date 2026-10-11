#!/usr/bin/env python3
"""Compare downloaded MaCRUISE outputs with saved PyMaCRUISE products.

This is a held-fixed, same-voxel-grid comparison. Source level-set NIfTIs
preserve voxel arrays but may have lost the original oblique affine; level-set
arrays are therefore compared by index. Source VTK points are voxel indices
and are transformed with the source Raw/T1 affine before RAS-mm surface
distances are measured.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.surfaces import read_legacy_vtk
from pymacruise.validation import (
    dice_coefficient,
    label_dice_scores,
    mesh_surface_distance_metrics,
    voxel_error_metrics,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path("/Users/landmaba/dev/mipav/sample-macruise/output")
DEFAULT_PYTHON = ROOT / "build/source-output-direct-comparison-2026-10-11/python"
DEFAULT_CORE = ROOT / "build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-2026-10-11"
DEFAULT_PREPARATION = ROOT / "build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-2026-10-11"
DEFAULT_REPORT = ROOT / "build/source-output-direct-comparison-2026-10-11/comparison.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def image_record(path: Path) -> dict[str, object]:
    image = nib.load(str(path))
    return {
        "path": str(path),
        "sha256": sha256(path),
        "shape": list(image.shape),
        "axis_codes": list(nib.aff2axcodes(image.affine)),
        "affine": image.affine.tolist(),
    }


def require_same_grid(reference: nib.spatialimages.SpatialImage, candidate: nib.spatialimages.SpatialImage, label: str) -> None:
    if reference.shape != candidate.shape or not np.allclose(reference.affine, candidate.affine, rtol=0, atol=1e-4):
        raise ValueError(f"{label} does not share the source T1 grid")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-output", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--python-output", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--python-core", type=Path, default=DEFAULT_CORE)
    parser.add_argument("--python-preparation", type=Path, default=DEFAULT_PREPARATION)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    source = args.source_output.resolve()
    python_output = args.python_output.resolve()
    core = args.python_core.resolve()
    preparation = args.python_preparation.resolve()
    source_t1_path = source / "Raw/T1.nii.gz"
    source_seg_path = source / "Raw/Seg.nii.gz"
    source_t1 = nib.load(str(source_t1_path))
    source_seg = nib.load(str(source_seg_path))
    python_input_t1 = nib.load(str(source.parent / "input/T1.nii.gz"))
    python_input_seg = nib.load(str(source.parent / "input/orig_target_seg.nii.gz"))
    require_same_grid(source_t1, source_seg, "Source Raw/Seg")
    if not np.array_equal(np.asanyarray(source_t1.dataobj), np.asanyarray(python_input_t1.dataobj)):
        raise ValueError("Source Raw/T1 does not match sample input T1 voxel-for-voxel")
    if not np.array_equal(np.asanyarray(source_seg.dataobj), np.asanyarray(python_input_seg.dataobj)):
        raise ValueError("Source Raw/Seg does not match sample input MaSeg voxel-for-voxel")

    report: dict[str, object] = {
        "scope": "same-grid held-fixed human source-output comparison; not whole-workflow source equivalence",
        "source_output": str(source),
        "input_identity": {
            "t1_voxel_identical": True,
            "maseg_voxel_identical": True,
            "t1": image_record(source_t1_path),
            "maseg": image_record(source_seg_path),
        },
        "settings_and_logs": {
            "source_run_logs_found": False,
            "note": "The output tree supplied for this comparison contains no run log or machine-readable run manifest; source runner parameters and exact build are therefore unverified.",
        },
        "tissue_masks": {},
        "levelsets": {},
        "segmentation_refinement": {},
        "surface_meshes": {},
    }

    source_mask_names = {
        "Atlas_GM": "atlas_gm.nii.gz",
        "Atlas_WM": "atlas_wm.nii.gz",
        "Atlas_Mask": "atlas_mask.nii.gz",
        "Atlas_BrainTissueMask": "atlas_brain_tissue_mask.nii.gz",
    }
    tissue_results = report["tissue_masks"]
    assert isinstance(tissue_results, dict)
    for source_name, python_name in source_mask_names.items():
        source_path = source / "CorticalSeperation" / f"{source_name}.nii.gz"
        python_path = preparation / python_name
        source_image, python_image = nib.load(str(source_path)), nib.load(str(python_path))
        require_same_grid(source_t1, source_image, f"Source {source_name}")
        require_same_grid(source_t1, python_image, f"Python {python_name}")
        source_values = np.asanyarray(source_image.dataobj) > 0
        python_values = np.asanyarray(python_image.dataobj) > 0
        tissue_results[source_name] = {
            "python_file": str(python_path),
            "source_positive_voxels": int(source_values.sum()),
            "python_positive_voxels": int(python_values.sum()),
            "binary_dice": dice_coefficient(source_values, python_values),
            "differing_voxels": int(np.count_nonzero(source_values != python_values)),
        }

    levels = report["levelsets"]
    assert isinstance(levels, dict)
    for name in ("inner", "central", "outer"):
        source_path = source / "Surfaces" / f"target_image_GMimg_{name}Vol.nii.gz"
        python_path = core / f"{name}.nii.gz"
        source_image, python_image = nib.load(str(source_path)), nib.load(str(python_path))
        if source_image.shape != python_image.shape or source_image.shape != source_t1.shape:
            raise ValueError(f"{name} level-set array shapes do not match the source voxel grid")
        source_values = np.asarray(source_image.dataobj, dtype=np.float32)
        python_values = np.asarray(python_image.dataobj, dtype=np.float32)
        levels[name] = {
            "source_file": str(source_path),
            "python_file": str(python_path),
            "source_affine_matches_raw_t1": bool(np.allclose(source_image.affine, source_t1.affine)),
            "python_affine_matches_raw_t1": bool(np.allclose(python_image.affine, source_t1.affine)),
            "voxel_error": voxel_error_metrics(source_values, python_values),
            "nonpositive_mask_dice": dice_coefficient(source_values <= 0, python_values <= 0),
            "sign_flip_voxels": int(np.count_nonzero((source_values <= 0) != (python_values <= 0))),
        }

    source_refine_path = source / "SegRefine/SEG_refine.nii.gz"
    python_refine_path = python_output / "SEG_refine.nii.gz"
    source_refined, python_refined = nib.load(str(source_refine_path)), nib.load(str(python_refine_path))
    require_same_grid(source_t1, source_refined, "Source SEG_refine")
    require_same_grid(source_t1, python_refined, "Python SEG_refine")
    source_labels = np.asarray(source_refined.dataobj, dtype=np.int32)
    python_labels = np.asarray(python_refined.dataobj, dtype=np.int32)
    per_label = label_dice_scores(source_labels, python_labels)
    dice_values = np.asarray(list(per_label.values()), dtype=np.float64)
    report["segmentation_refinement"] = {
        "source_file": str(source_refine_path),
        "python_file": str(python_refine_path),
        "differing_voxels": int(np.count_nonzero(source_labels != python_labels)),
        "nonbackground_label_count": len(per_label),
        "minimum_dice": float(dice_values.min()),
        "mean_dice": float(dice_values.mean()),
        "median_dice": float(np.median(dice_values)),
        "labels_below_0_97": int(np.count_nonzero(dice_values < 0.97)),
        "per_label_dice": {str(label): score for label, score in per_label.items()},
        "labels_only_in_source": [int(label) for label in per_label if np.any(source_labels == label) and not np.any(python_labels == label)],
        "labels_only_in_python": [int(label) for label in per_label if np.any(python_labels == label) and not np.any(source_labels == label)],
    }

    meshes = report["surface_meshes"]
    assert isinstance(meshes, dict)
    source_affine = source_t1.affine
    for name in ("inner", "central", "outer"):
        source_path = source / "Surfaces" / f"target_image_GMimg_{name}Surf.vtk"
        python_path = python_output / "surfaces" / f"{name}_surface.vtk"
        source_mesh = read_legacy_vtk(source_path, input_frame="RAS")
        python_mesh = read_legacy_vtk(python_path, input_frame="RAS")
        source_points = (np.c_[source_mesh.vertices, np.ones(len(source_mesh.vertices))] @ source_affine.T)[:, :3]
        distances = mesh_surface_distance_metrics(
            source_points, source_mesh.faces, python_mesh.vertices, python_mesh.faces
        )
        meshes[name] = {
            "source_file": str(source_path),
            "python_file": str(python_path),
            "source_point_units_before_transform": "voxel indices",
            "source_transform": source_affine.tolist(),
            "source_topology": {"vertices": len(source_mesh.vertices), "faces": len(source_mesh.faces)},
            "python_topology": {"vertices": len(python_mesh.vertices), "faces": len(python_mesh.faces)},
            "symmetric_vertex_to_triangle_mm": distances,
            "passes_user_surface_limits": bool(
                distances["mean_absolute_distance_mm"] <= 0.01 and distances["hausdorff_mm"] <= 1.0
            ),
        }

    source_thickness_path = source / "Thickness/target_image_GMimg_innerVol_thick.nii.gz"
    python_thickness_path = python_output / "surfaces/thickness_levelset.nii.gz"
    source_thickness, python_thickness = nib.load(str(source_thickness_path)), nib.load(str(python_thickness_path))
    if source_thickness.shape != python_thickness.shape or source_thickness.shape != source_t1.shape:
        raise ValueError("Thickness arrays do not share the source voxel grid")
    source_values = np.asarray(source_thickness.dataobj, dtype=np.float32)
    python_values = np.asarray(python_thickness.dataobj, dtype=np.float32)
    common_positive = (source_values > 0) & (python_values > 0)
    delta = np.abs(source_values[common_positive] - python_values[common_positive])
    report["thickness_map"] = {
        "source_file": str(source_thickness_path),
        "python_file": str(python_thickness_path),
        "source_affine_matches_raw_t1": bool(np.allclose(source_thickness.affine, source_t1.affine)),
        "python_affine_matches_raw_t1": bool(np.allclose(python_thickness.affine, source_t1.affine)),
        "source_positive_voxels": int((source_values > 0).sum()),
        "python_positive_voxels": int((python_values > 0).sum()),
        "positive_mask_dice": dice_coefficient(source_values > 0, python_values > 0),
        "common_positive_voxels": int(common_positive.sum()),
        "common_positive_mae_mm": float(delta.mean()) if delta.size else None,
        "common_positive_rmse_mm": float(np.sqrt(np.mean(np.square(delta)))) if delta.size else None,
        "common_positive_max_abs_mm": float(delta.max()) if delta.size else None,
    }

    destination = args.report.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(destination), "surface_gates_pass": all(v["passes_user_surface_limits"] for v in meshes.values()), "minimum_refinement_label_dice": report["segmentation_refinement"]["minimum_dice"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
