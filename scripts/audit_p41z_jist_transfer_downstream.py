#!/usr/bin/env python3
"""Run selected plus surface stages with Java-compatible cropped Dura maps.

This is a synthetic compatibility probe. It holds the P41z crop-ready
segmentation and atlas maps fixed, replaces only TOADS Dura maps with the
selected Java transfer-plus-crop arrays, and runs topology through surfaces.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cropping import CropParameters
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_post_tgdm_stages
from pymacruise.plus_topology_composition import run_plus_topology_and_tgdm
from pymacruise.surfaces import SurfaceMesh, write_legacy_vtk
from pymacruise.tgdm import TGDMMaCRUISEImageResult
from pymacruise.validation import mesh_topology_metrics

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "build/p41z-case-normalization/iterations-0008-to-0011-surfaces"
JAVA_CROP = ROOT / "build/p41z-case-normalization/iteration-0011-java-crop-transfer/java-output"
MAPS = ("csf", "gm", "wm_fill", "wm_mask")
INPUTS = (
    "segmentation", "toads_wm_mask", "multi_atlas_gm", "multi_atlas_wm",
    "multi_atlas_mask", "toads_gm", "toads_wm_fill", "toads_wm_mask",
)


def _load(path: Path) -> nib.Nifti1Image:
    image = nib.load(str(path))
    if not isinstance(image, nib.Nifti1Image) or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"expected RAS NIfTI-1 image: {path}")
    return image


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/p41z-jist-transfer-downstream")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    plus_root = CASE / "plus-surfaces"
    crop_attempt = plus_root / "stages/transfer_crop/attempt-01"
    crop_record = json.loads((crop_attempt / "crop.json").read_text())
    crop = CropParameters(
        original_shape=tuple(crop_record["original_shape"]),
        bounds=tuple(tuple(pair) for pair in crop_record["bounds"]),
        background=float(crop_record["background"]),
        original_affine=np.asarray(crop_record["original_affine"], dtype=np.float64),
    )
    images = {name: _load(crop_attempt / f"{name}.nii.gz") for name in INPUTS}
    shape = images["segmentation"].shape
    if any(image.shape != shape or not np.allclose(image.affine, images["segmentation"].affine, atol=1e-5, rtol=0) for image in images.values()):
        raise ValueError("saved crop-ready maps do not share one RAS grid")

    raw_shape = tuple(int(v) for v in (JAVA_CROP / "bounds.txt").read_text().strip().split(";")[-1].split(","))
    if raw_shape != shape:
        raise ValueError(f"Java crop shape {raw_shape} differs from prepared crop {shape}")
    java_maps: dict[str, nib.Nifti1Image] = {}
    for name in MAPS:
        values = np.fromfile(JAVA_CROP / f"{name}.f32be", dtype=">f4").astype(np.float32).reshape(raw_shape, order="C")
        java_maps[name] = nib.Nifti1Image(values, images["segmentation"].affine)
        nib.save(java_maps[name], output / f"jist_cropped_{name}.nii.gz")

    # Dura CSF is carried for provenance but the selected downstream plus
    # composition routes GM, WM fill, and discrete WM mask only.
    source_inputs = (
        images["segmentation"], java_maps["wm_mask"], images["multi_atlas_gm"],
        images["multi_atlas_wm"], images["multi_atlas_mask"], java_maps["gm"],
        java_maps["wm_fill"], java_maps["wm_mask"],
    )
    composed = run_plus_topology_and_tgdm(*source_inputs, options=PlusPipelineOptions())
    tgdm = composed.tgdm
    post = run_plus_post_tgdm_stages(tgdm, crop)
    level_sets = {
        "inner": post.level_sets.inner,
        "central": post.level_sets.central,
        "outer": post.level_sets.outer,
    }
    for name, image in level_sets.items():
        nib.save(image, output / f"{name}_level_set.nii.gz")
    meshes: dict[str, SurfaceMesh] = {
        "inner": post.surfaces.inner,
        "central": post.surfaces.central,
        "outer": post.surfaces.outer,
    }
    for name, mesh in meshes.items():
        write_legacy_vtk(mesh, output / f"{name}_surface.vtk")

    aseg = _load(CASE / "maseg-cortical-separation/Atlas_GM.nii.gz")
    labels = np.asarray(_load(CASE / "maseg-cortical-separation/Atlas_BrainTissueMask.nii.gz").dataobj)
    inner = np.asarray(level_sets["inner"].dataobj) <= 0
    outer = np.asarray(level_sets["outer"].dataobj) <= 0
    if labels.shape != inner.shape:
        raise ValueError("full-grid MaSeg tissue labels and restored level sets differ")
    wm = labels == 44
    gm = labels == 100
    csf = labels == 5
    support = {
        "wm_inside_inner": float(np.count_nonzero(wm & inner) / max(1, np.count_nonzero(wm))),
        "gm_between_inner_outer": float(np.count_nonzero(gm & ~inner & outer) / max(1, np.count_nonzero(gm))),
        "csf_outside_outer": float(1 - np.count_nonzero(csf & outer) / max(1, np.count_nonzero(csf))),
    }
    topology = {
        name: mesh_topology_metrics(mesh.vertices, mesh.faces)
        for name, mesh in meshes.items()
    }
    physical_quality = json.loads((CASE / "ribbon-quality.json").read_text())
    physical_support = physical_quality["label_surface_consistency"]["support_fractions"]
    support_key_map = {
        "wm_inside_inner": "wm_inside_inner",
        "gm_between_inner_outer": "gm_between_inner_and_outer",
        "csf_outside_outer": "csf_outside_outer",
    }
    physical_surfaces = plus_root / "stages/surfaces_thickness/attempt-01"
    baseline_changes = {}
    for name, image in level_sets.items():
        baseline = np.asarray(_load(physical_surfaces / f"{name}_level_set.nii.gz").dataobj)
        values = np.asarray(image.dataobj)
        baseline_changes[name] = {
            "sign_differences": int(np.count_nonzero((values <= 0) != (baseline <= 0))),
            "mean_absolute_levelset_difference": float(np.mean(np.abs(values.astype(np.float32) - baseline.astype(np.float32)), dtype=np.float64)),
        }

    report = {
        "evidence_category": "synthetic downstream compatibility experiment; selected Java crop inputs, Python topology-through-surface stages",
        "source_transfer_behavior": "JIST-preserved Dura extents and selected Java crop",
        "p41z_checkpoint": "LesionTOADS iteration 11; nonterminal upstream checkpoint",
        "source_crop_grid_shape": list(shape),
        "source_crop_bounds_inclusive": [list(pair) for pair in crop.bounds],
        "input_sha256": {
            "runner": _sha256(Path(__file__).resolve()),
            "crop_parameters": _sha256(crop_attempt / "crop.json"),
            "source_compatible_java_crop_maps": {
                name: _sha256(JAVA_CROP / f"{name}.f32be") for name in MAPS
            },
            "fixed_crop_ready_inputs": {
                name: _sha256(crop_attempt / f"{name}.nii.gz") for name in INPUTS
            },
            "physical_route_ribbon_reference": _sha256(CASE / "ribbon-quality.json"),
        },
        "compared_with": "saved physical-affine transfer plus crop surface run",
        "physical_route_quality_reference": physical_support,
        "jist_crop_route_label_support": support,
        "label_support_delta_vs_physical_percentage_points": {
            name: 100 * (support[name] - physical_support[support_key_map[name]])
            for name in support
        },
        "levelset_changes_vs_physical_route": baseline_changes,
        "surface_topology": topology,
        "thickness": {
            "surface_mean_mm": float(post.thickness.surface.mean_mm),
            "surface_median_mm": float(np.median(post.thickness.surface.inner_surface.point_data["mean_thickness_mm"])),
            "surface_p05_mm": float(np.percentile(post.thickness.surface.inner_surface.point_data["mean_thickness_mm"], 5)),
            "surface_p95_mm": float(np.percentile(post.thickness.surface.inner_surface.point_data["mean_thickness_mm"], 95)),
            "levelset_mean_mm": float(post.thickness.volume.mean_mm),
        },
        "passes_97_percent_tissue_support": bool(min(support.values()) >= 0.97),
        "source_scope_limit": "This is not a Java surface evolution comparison, a full PLUS run, or real-MRI validation.",
    }
    (output / "downstream-audit.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"report={output / 'downstream-audit.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
