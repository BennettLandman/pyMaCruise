#!/usr/bin/env python3
"""Compare an independent selected-source PLUS run with a Python full run.

The external runner's output layout is deliberately supplied as a JSON product
map because the pinned runtime's fresh PLUS export layout has not been tested.
This tool compares products; it cannot establish whole-port acceptance when
source identity, preprocessing boundaries, or stage outputs are incomplete.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from scipy.io import loadmat

from pymacruise.surfaces import read_legacy_vtk
from pymacruise.validation import (
    dice_coefficient, label_dice_scores, mesh_surface_distance_metrics,
    voxel_error_metrics,
)


VOLUMES = ("segmentation", "inner_level_set", "central_level_set", "outer_level_set", "thickness_volume")
MESHES = ("inner_surface", "central_surface", "outer_surface")
PRODUCTS = (*VOLUMES, *MESHES, "surface_thickness_mat")
REQUIRED_BASE_STAGES = {
    "preprocessing", "crop", "initial_alignment", "registered_recrop",
    "lesion_workflow", "dura_outputs", "native_dura_transfer",
    "maseg_cortical_separation", "plus_surface_workflow", "maseg_refinement",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _regular(path: Path, description: str) -> Path:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Missing regular {description}: {path}")
    return path.resolve()


def _records(value):
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
            yield value
        else:
            for child in value.values():
                yield from _records(child)
    elif isinstance(value, list):
        for child in value:
            yield from _records(child)


def verify_hash_records(document: dict, context: str) -> int:
    count = 0
    for record in _records(document):
        path = _regular(Path(record["path"]), f"{context} output")
        if sha256(path) != record["sha256"]:
            raise ValueError(f"{context} output hash mismatch: {path}")
        count += 1
    assets = document.get("inputs", {}).get("atlasem_assets", {})
    if isinstance(assets, dict):
        for path_text, expected in assets.items():
            path = _regular(Path(path_text), f"{context} AtlasEM payload")
            if sha256(path) != expected:
                raise ValueError(f"{context} AtlasEM payload hash mismatch: {path}")
            count += 1
    return count


def _image(path: Path):
    path = _regular(path, "NIfTI product")
    image = nib.load(str(path))
    data = np.asanyarray(image.dataobj)
    if (not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) or len(image.shape) != 3
            or np.iscomplexobj(data) or not np.issubdtype(data.dtype, np.number)
            or not np.isfinite(data).all() or not np.isfinite(image.affine).all()):
        raise ValueError(f"Expected finite 3D NIfTI: {path}")
    spatial_units, _ = image.header.get_xyzt_units()
    if spatial_units != "mm" or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"NIfTI must use RAS axes and millimetre units: {path}")
    return image, data


def _distribution(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=np.float64).ravel()
    if not values.size or not np.isfinite(values).all():
        raise ValueError("Thickness values must be nonempty and finite")
    return {"count": int(values.size), "mean_mm": float(values.mean()),
            "std_mm": float(values.std()), "min_mm": float(values.min()),
            "median_mm": float(np.median(values)), "max_mm": float(values.max())}


def _python_products(manifest: dict) -> dict[str, Path]:
    stages = manifest["stages"]
    surface_manifest_path = Path(stages["plus_surface_workflow"]["manifest"]["path"])
    surface_manifest = json.loads(_regular(surface_manifest_path, "surface child manifest").read_text())
    if surface_manifest.get("status") != "complete":
        raise ValueError("PLUS surface/thickness child manifest is not complete")
    outputs = surface_manifest["stages"]["surfaces_thickness"]["outputs"]
    refinement_manifest_path = Path(stages["maseg_refinement"]["manifest"]["path"])
    refinement = json.loads(_regular(refinement_manifest_path, "refinement child manifest").read_text())
    if refinement.get("status") != "complete":
        raise ValueError("Final MaSeg refinement child manifest is not complete")
    refined = refinement["outputs"]["SEG_refine.nii.gz"]["path"]
    return {
        "segmentation": Path(refined),
        "inner_level_set": Path(outputs["inner_level_set"]["path"]),
        "central_level_set": Path(outputs["central_level_set"]["path"]),
        "outer_level_set": Path(outputs["outer_level_set"]["path"]),
        "thickness_volume": Path(outputs["thickness_volume"]["path"]),
        "inner_surface": Path(outputs["inner_surface"]["path"]),
        "central_surface": Path(outputs["central_surface"]["path"]),
        "outer_surface": Path(outputs["outer_surface"]["path"]),
        "surface_thickness_mat": Path(outputs["thickness_metrics"]["path"]),
    }


def compare(args) -> dict:
    atlas_policy = getattr(args, "atlas_policy", "strict")
    atlas_differences = []
    def atlas_check(matches: bool, message: str) -> None:
        if matches:
            return
        if atlas_policy == "strict":
            raise ValueError(message)
        atlas_differences.append(message)
    run_path = _regular(args.python_run / "plus-full-workflow.json", "Python PLUS manifest")
    manifest = json.loads(run_path.read_text())
    if manifest.get("status") != "complete" or not manifest.get("terminal"):
        raise ValueError("Python PLUS run must be complete and terminal")
    atlasem_enabled = bool(manifest.get("settings", {}).get("atlasem_dura_enabled"))
    if not atlasem_enabled:
        raise ValueError("Selected complete PLUS comparison requires the AtlasEM/Dura stage")
    expected = REQUIRED_BASE_STAGES | ({"atlasem_dura"} if atlasem_enabled else set())
    expected_order = ["preprocessing"]
    if atlasem_enabled:
        expected_order.append("atlasem_dura")
    expected_order.extend(("crop", "initial_alignment", "registered_recrop", "lesion_workflow",
                           "dura_outputs", "native_dura_transfer", "maseg_cortical_separation",
                           "plus_surface_workflow", "maseg_refinement"))
    order = manifest.get("stage_order", [])
    if order != expected_order:
        raise ValueError("Python PLUS stage order differs from the selected workflow contract")
    if set(order) != expected or set(manifest.get("stages", {})) != expected:
        raise ValueError("Python PLUS manifest has missing, extra, or ambiguous active stages")
    verified_records = verify_hash_records(manifest, "Python PLUS parent")
    for stage, record in manifest["stages"].items():
        if record.get("status") != "complete":
            raise ValueError(f"Python PLUS stage is incomplete: {stage}")
        child_record = record.get("manifest", {})
        child_path = _regular(Path(child_record.get("path", "")), f"{stage} child manifest")
        if sha256(child_path) != child_record.get("sha256"):
            raise ValueError(f"Python PLUS child manifest hash mismatch: {stage}")
        child = json.loads(child_path.read_text())
        if child.get("status") != "complete":
            raise ValueError(f"Python PLUS child manifest is incomplete: {stage}")
        verified_records += verify_hash_records(child, stage)
    surface_manifest = json.loads(Path(manifest["stages"]["plus_surface_workflow"]["manifest"]["path"]).read_text())
    for stage, record in surface_manifest["stages"].items():
        if record.get("status") != "complete":
            raise ValueError(f"Python PLUS surface substage is incomplete: {stage}")
        verified_records += verify_hash_records(record, f"surface/{stage}")

    info_path = _regular(args.reference_info, "external reference run record")
    info = json.loads(info_path.read_text())
    if info.get("status") != "complete" or info.get("exit_code") != 0:
        raise ValueError("External reference record must declare a complete run and zero exit status")
    if not info.get("runtime") or not info.get("source_revision") or not info.get("command"):
        raise ValueError("External reference record must identify runtime, source revision, and launch command")
    log_path = _regular(args.reference_log, "external launch log")
    if log_path.stat().st_size == 0:
        raise ValueError("External launch log is empty")
    input_paths = {"t1": args.t1, "flair": args.flair, "brain_mask": args.brain_mask, "maseg": args.maseg}
    input_hashes = {name: sha256(_regular(path, f"supplied {name}")) for name, path in input_paths.items()}
    if any(info.get("input_sha256", {}).get(name) != digest
           for name, digest in input_hashes.items()):
        raise ValueError("External run input hashes do not match the declared independent input packet")
    for name, expected_hash in input_hashes.items():
        record = manifest.get("inputs", {}).get(name)
        if not isinstance(record, dict) or record.get("sha256") != expected_hash:
            raise ValueError(f"Python PLUS run did not use the same {name} input")
    if info.get("atlasem_enabled") is not True:
        raise ValueError("External and Python AtlasEM stage selection differs")
    if info.get("products_scope") != "independently_generated_selected_plus_products":
        raise ValueError("External products must be declared independently generated selected PLUS outputs")
    source_mesh_frame = getattr(args, "source_mesh_frame", "ras-mm")
    expected_frame = "RAS+ mm" if source_mesh_frame == "ras-mm" else "T1 voxel indices"
    if not info.get("preprocessing_boundary") or info.get("surface_frame") != expected_frame:
        raise ValueError(f"External run record must declare its preprocessing boundary and mesh frame {expected_frame}")
    required_labels = info.get("required_segmentation_labels")
    if not isinstance(required_labels, list) or not required_labels or any(
        isinstance(label, bool) or not isinstance(label, int) or label == 0 for label in required_labels
    ):
        raise ValueError("External run record must list the required nonbackground segmentation labels")
    attempts = info.get("active_product_attempts")
    if (not isinstance(attempts, dict) or set(attempts) != set(PRODUCTS)
            or any(not isinstance(attempt, str) or not attempt.strip() for attempt in attempts.values())):
        raise ValueError("External run record must identify the active attempt/revision for every product")
    for name, record in manifest.get("inputs", {}).items():
        if name in input_hashes:
            continue  # Four scientific starting inputs were checked strictly above.
        if info.get("input_sha256", {}).get(name) != record.get("sha256"):
            atlas_check(False, f"External and Python atlas input identity differs for {name}")
    atlasem_assets: dict[str, str] | None = None
    if atlasem_enabled:
        atlasem_path = Path(manifest["stages"]["atlasem_dura"]["manifest"]["path"])
        atlasem_manifest = json.loads(_regular(atlasem_path, "AtlasEM/Dura child manifest").read_text())
        atlasem_descriptor = atlasem_manifest["inputs"]["atlasem_descriptor"]["sha256"]
        raw_assets = atlasem_manifest["inputs"].get("atlasem_assets", {})
        atlasem_assets = {}
        for name, digest in raw_assets.items():
            key = Path(name).name
            if key in atlasem_assets:
                raise ValueError("AtlasEM payload filenames are ambiguous; preserve relative asset names in run records")
            atlasem_assets[key] = digest
        external_assets = info.get("atlasem_asset_sha256")
        atlas_check(not (info.get("atlasem_descriptor_sha256") != atlasem_descriptor
                or not isinstance(external_assets, dict)
                or external_assets != atlasem_assets),
                "External AtlasEM descriptor or payload hashes differ from the Python run")
    alignment_path = Path(manifest["stages"]["initial_alignment"]["manifest"]["path"])
    alignment_manifest = json.loads(_regular(alignment_path, "LesionTOADS alignment manifest").read_text())
    raw_shape_assets = alignment_manifest["inputs"].get("shape_atlases", {})
    lesion_atlas_assets: dict[str, str] = {}
    for record in raw_shape_assets.values():
        key = Path(record["path"]).name
        if key in lesion_atlas_assets:
            raise ValueError("LesionTOADS atlas payload filenames are ambiguous; retain relative asset names")
        lesion_atlas_assets[key] = record["sha256"]
    atlas_check(bool(lesion_atlas_assets) and info.get("lesion_atlas_asset_sha256") == lesion_atlas_assets,
                "External LesionTOADS descriptor payload hashes differ from the Python run")
    source_map = json.loads(_regular(args.product_map, "external product map").read_text())
    if set(source_map) != set(PRODUCTS):
        raise ValueError(f"External product map must contain exactly: {', '.join(PRODUCTS)}")
    external_root = args.external_root.resolve()
    external = {}
    for name, value in source_map.items():
        path = Path(value)
        path = path.resolve() if path.is_absolute() else (external_root / path).resolve()
        if path != external_root and external_root not in path.parents:
            raise ValueError(f"External product escapes declared output root: {name}")
        external[name] = _regular(path, f"external {name}")
    if len(set(external.values())) != len(external):
        raise ValueError("External product map aliases multiple required products to the same file")
    recorded_product_hashes = info.get("product_sha256")
    actual_product_hashes = {name: sha256(path) for name, path in external.items()}
    if recorded_product_hashes != actual_product_hashes:
        raise ValueError("External run record product hashes do not match the active product map")
    python = _python_products(manifest)
    if not set(python) == set(PRODUCTS):
        raise ValueError("Python PLUS products are incomplete")

    t1_img, _ = _image(args.t1)
    report = {
        "scope": "independent selected PLUS product comparison",
        "overall_port_acceptance": "unverified",
        "atlas_identity_policy": atlas_policy,
        "atlas_identity_differences": atlas_differences,
        "source_mesh_frame": source_mesh_frame,
        "source_revision": info["source_revision"], "runtime": info["runtime"],
        "external_run_record": str(info_path), "external_run_record_sha256": sha256(info_path),
        "external_log": str(log_path), "external_log_sha256": sha256(log_path),
        "external_root": str(external_root), "python_manifest": str(run_path),
        "python_manifest_sha256": sha256(run_path), "verified_python_records": verified_records,
        "input_sha256": input_hashes, "atlasem_enabled": atlasem_enabled,
        "atlasem_asset_sha256": atlasem_assets,
        "lesion_atlas_asset_sha256": lesion_atlas_assets,
        "preprocessing_boundary": info["preprocessing_boundary"],
        "required_segmentation_labels": required_labels,
        "active_external_product_attempts": attempts,
        "products": {}, "limitations": [
            "The selected AtlasEM atlas payload identity and source/runtime settings must be independently audited.",
            "No universal scalar or thickness acceptance threshold is established.",
            "Surface thickness has no established vertex correspondence across independently generated meshes.",
            "This product adapter does not validate every upstream intermediate or transform.",
        ],
    }
    images = {}
    for name in VOLUMES:
        ext_img, ext_values = _image(external[name])
        py_img, py_values = _image(python[name])
        if any(img.shape != t1_img.shape for img in (ext_img, py_img)):
            raise ValueError(f"{name} must be on the full T1 input grid")
        if not np.allclose(ext_img.affine, py_img.affine, atol=1e-4, rtol=0):
            raise ValueError(f"{name} external and Python physical grids differ")
        if not np.allclose(ext_img.affine, t1_img.affine, atol=1e-4, rtol=0):
            raise ValueError(f"{name} does not use the declared native T1 RAS frame")
        images[name] = (ext_values, py_values)
        entry = {"external": str(external[name]), "external_sha256": sha256(external[name]),
                 "python": str(python[name]), "python_sha256": sha256(python[name]),
                 "affine_ras_mm": ext_img.affine.tolist()}
        if name == "segmentation":
            if not np.array_equal(ext_values, np.rint(ext_values)) or not np.array_equal(py_values, np.rint(py_values)):
                raise ValueError("Final segmentation labels must be integers")
            scores = label_dice_scores(ext_values, py_values, labels=tuple(required_labels))
            if not scores:
                raise ValueError("Final segmentation contains no required nonbackground labels")
            entry.update(per_label_dice={str(k): float(v) for k, v in scores.items()},
                         minimum_dice=float(min(scores.values())),
                         labels_below_threshold=[int(k) for k, v in scores.items() if v < .97],
                         differing_voxels=int(np.count_nonzero(ext_values != py_values)),
                         passed=all(v >= .97 for v in scores.values()))
        else:
            entry["voxel_error"] = voxel_error_metrics(ext_values, py_values)
            entry["negative_region_dice"] = dice_coefficient(ext_values <= 0, py_values <= 0)
            if name == "thickness_volume":
                for role, values in (("external", ext_values), ("python", py_values)):
                    positive = np.asarray(values)[np.asarray(values) > 0]
                    entry[f"{role}_positive_thickness_distribution"] = (
                        _distribution(positive) if positive.size else None
                    )
        report["products"][name] = entry

    for name in MESHES:
        ext = read_legacy_vtk(external[name], input_frame="RAS")
        py = read_legacy_vtk(python[name], input_frame="RAS")
        ext_vertices = ext.vertices
        if source_mesh_frame == "voxel":
            ext_vertices = nib.affines.apply_affine(t1_img.affine, ext_vertices)
        distances = mesh_surface_distance_metrics(ext_vertices, ext.faces, py.vertices, py.faces)
        report["products"][name] = {
            "external": str(external[name]), "external_sha256": sha256(external[name]),
            "python": str(python[name]), "python_sha256": sha256(python[name]),
            **distances,
            "passed": distances["mean_absolute_distance_mm"] <= .01 and distances["hausdorff_mm"] <= 1.,
        }
    ext_mat = loadmat(external["surface_thickness_mat"])
    py_metrics = json.loads(python["surface_thickness_mat"].read_text())
    if "whole_thickness" not in ext_mat or "surface" not in py_metrics:
        raise ValueError("Thickness summary is missing required whole_thickness/surface data")
    ext_thickness = np.asarray(ext_mat["whole_thickness"]).ravel()
    # Current Python reports aggregate thickness metrics, not per-vertex samples.
    report["products"]["surface_thickness"] = {
        "external": str(external["surface_thickness_mat"]),
        "external_sha256": sha256(external["surface_thickness_mat"]),
        "python": str(python["surface_thickness_mat"]),
        "python_sha256": sha256(python["surface_thickness_mat"]),
        "external_distribution": _distribution(ext_thickness),
        "python_summary": py_metrics["surface"],
        "summary_absolute_errors": {
            "mean_mm": abs(float(py_metrics["surface"]["mean_mm"]) - float(np.mean(ext_thickness))),
            "standard_deviation_mm": abs(float(py_metrics["surface"]["standard_deviation_mm"]) - float(np.std(ext_thickness))),
            "sample_count": abs(int(py_metrics["surface"]["nonzero_vertex_count"]) - int(ext_thickness.size)),
        },
        "vertexwise_comparison": "unavailable: the current Python thickness manifest stores aggregates only",
    }
    report["segmentation_and_surface_gates_passed"] = (
        report["products"]["segmentation"]["passed"]
        and all(report["products"][name]["passed"] for name in MESHES)
    )
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("external-root", "python-run", "t1", "flair", "brain-mask", "maseg",
                 "reference-info", "reference-log", "product-map", "report"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--atlas-policy", choices=("strict", "report-differences"), default="strict",
                        help="strict rejects different atlas bytes; report-differences permits exploratory product metrics and explicitly leaves atlas identity unresolved")
    parser.add_argument("--source-mesh-frame", choices=("voxel", "ras-mm"), default="ras-mm",
                        help="explicit source mesh coordinates; voxel applies the declared T1 affine without fitting")
    args = parser.parse_args(argv)
    try:
        report = compare(args)
        status = 0 if report["segmentation_and_surface_gates_passed"] else 1
    except (ValueError, KeyError, OSError, TypeError) as exc:
        report = {"overall_port_acceptance": "unverified", "comparison_error": str(exc)}
        status = 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Report: {args.report}; exit code: {status}", flush=True)
    return status


if __name__ == "__main__":
    sys.exit(main())
