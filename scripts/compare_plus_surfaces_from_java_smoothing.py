#!/usr/bin/env python3
"""Diagnostic: hold the Java smoothing output fixed and compare later Python TGDM surfaces."""
from __future__ import annotations

import argparse
import gzip
import json
import struct
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.affines import apply_affine
from skimage.measure import marching_cubes

from pymacruise.surfaces import extract_levelset_surface
from pymacruise.tgdm import (
    TGDMBaselineParameters,
    _INNER_FROM_CENTRAL_PROFILE_SETTINGS,
    tgdm_central_pressure_field,
    tgdm_evolve,
    tgdm_inner_pressure_field,
    tgdm_outer_pressure_field,
)
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics


def _read_java_volume(root: Path, name: str, shape: tuple[int, int, int]) -> np.ndarray:
    with gzip.open(root / "java" / f"{name}.raw.gz", "rb") as stream:
        return np.frombuffer(stream.read(), dtype=">f4").astype(np.float32).reshape(shape)


def _surface_metrics(
    java: np.ndarray, python: np.ndarray, affine: np.ndarray
) -> dict[str, float | int | bool]:
    java_surface = extract_levelset_surface(nib.Nifti1Image(java, affine))
    py_surface = extract_levelset_surface(nib.Nifti1Image(python, affine))
    java_vertices, java_faces, _, _ = marching_cubes(java, level=0.0)
    py_vertices, py_faces, _, _ = marching_cubes(python, level=0.0)
    metrics: dict[str, float | int | bool] = {
        "sign_differences": int(np.count_nonzero((java <= 0.0) != (python <= 0.0))),
        "dice_at_zero": dice_coefficient(java <= 0.0, python <= 0.0),
    }
    metrics.update(mesh_surface_distance_metrics(
        java_surface.vertices, java_surface.faces,
        py_surface.vertices, py_surface.faces,
    ))
    metrics["marching_cubes_diagnostic"] = mesh_surface_distance_metrics(
        apply_affine(affine, java_vertices), java_faces,
        apply_affine(affine, py_vertices), py_faces,
    )
    metrics["thresholds_pass"] = bool(
        metrics["dice_at_zero"] >= 0.97
        and metrics["mean_absolute_distance_mm"] <= 0.01
        and metrics["hausdorff_mm"] <= 1.0
    )
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    root = args.reference_dir.resolve()
    provenance = json.loads((root / "input-provenance.json").read_text())
    shape = tuple(provenance["shape_xyz"])
    voxel_count = int(np.prod(shape))
    with (root / "arrays.bin").open("rb") as stream:
        actual_shape = struct.unpack(">III", stream.read(12))
        if actual_shape != shape:
            raise SystemExit("held-fixed arrays do not match the provenance shape")
        stream.read(voxel_count * 4)  # initial membership
        fields = [
            np.frombuffer(stream.read(voxel_count * 4), dtype=">f4")
            .astype(np.float32).reshape(shape)
            for _ in range(4)
        ]
        gvf = np.frombuffer(stream.read(voxel_count * 12), dtype=">f4")
        gvf = gvf.astype(np.float32).reshape((*shape, 3))
        if stream.read(1):
            raise SystemExit("held-fixed array binary has unexpected trailing bytes")
    wm, gm, ace_gm, skeleton = fields
    java_smoothing = _read_java_volume(root, "smoothed_initial", shape)
    parameters = TGDMBaselineParameters()
    inner_force = tgdm_inner_pressure_field(
        wm, white_matter_iso=parameters.white_matter_iso
    )
    initial_inner = tgdm_evolve(
        java_smoothing,
        profile="inner",
        curvature_force=-parameters.inner_curvature,
        pressure_force=parameters.inner_pressure,
        external_force=parameters.inner_external,
        region_force=inner_force,
        outer_iterations=parameters.inner_iterations,
    )
    tissue_sign = np.load(root / "python-tissue_sign.npy").astype(np.uint8)
    central = tgdm_evolve(
        initial_inner.level_set,
        profile="central",
        curvature_force=-parameters.central_curvature,
        pressure_force=parameters.central_pressure,
        external_force=parameters.central_external,
        region_force=tgdm_central_pressure_field(tissue_sign),
        advection=gvf,
        inner_bound=initial_inner.level_set,
        outer_iterations=parameters.central_iterations,
    )
    central_shrink = tgdm_evolve(
        central.level_set,
        profile="inner",
        curvature_force=-parameters.central_curvature,
        pressure_force=-1.0,
        outer_iterations=parameters.central_shrink_iterations,
        settings=_INNER_FROM_CENTRAL_PROFILE_SETTINGS,
    )
    inner = tgdm_evolve(
        central_shrink.level_set,
        profile="inner",
        curvature_force=-parameters.inner_curvature,
        pressure_force=-parameters.inner_pressure,
        external_force=parameters.inner_external,
        region_force=-inner_force,
        outer_iterations=parameters.inner_iterations,
    )
    outer_force = tgdm_outer_pressure_field(
        central.level_set, gm, skeleton, gray_matter_iso=parameters.gray_matter_iso
    )
    outer = tgdm_evolve(
        central.level_set,
        profile="outer",
        curvature_force=-parameters.outer_curvature,
        pressure_force=parameters.outer_pressure,
        external_force=parameters.outer_external,
        region_force=outer_force,
        inner_bound=central.level_set,
        skeletal_bound=skeleton,
        outer_iterations=parameters.central_iterations,
    )
    affine = np.asarray(provenance["affine_ras"], dtype=np.float64)
    stages = {
        "initial_inner": initial_inner.level_set,
        "central": central.level_set,
        "inner": inner.level_set,
        "outer": outer.level_set,
    }
    results = {
        name: _surface_metrics(_read_java_volume(root, name, shape), field, affine)
        for name, field in stages.items()
    }
    report = {
        "diagnostic_only": True,
        "synthetic_only": provenance["synthetic_only"],
        "evidence_category": provenance["evidence_category"],
        "comparison": "selected Java TGDM outputs vs Python downstream evolution with Java smoothed_initial substituted",
        "scope": "tests stages after initial smoothing only; it is not an end-to-end Python-vs-Java result",
        "workflow_run_signature": provenance["workflow_run_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "surface_extractor": "source-matched JIST IsoSurfaceOnGrid connectivity-18 port",
        "surface_gate_note": "Primary distance metrics use production JIST-compatible surfaces; marching_cubes_diagnostic is informational only.",
        "affine_ras": affine.tolist(),
        "settings": parameters.__dict__,
        "outputs": results,
        "surface_thresholds_pass": all(bool(item["thresholds_pass"]) for item in results.values()),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    if not report["surface_thresholds_pass"]:
        raise SystemExit("one or more downstream held-fixed TGDM surfaces missed an agreed tolerance")


if __name__ == "__main__":
    main()
