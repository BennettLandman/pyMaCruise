#!/usr/bin/env python3
"""Compare a saved brain-resolution Dura stage against selected Java outputs.

Inputs are the exact four RAS NIfTIs passed to the Java Dura harness. This
allows source comparison against a current workflow state without changing or
overwriting the committed baseline fixtures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_dura_base_java import read_java_volume  # noqa: E402
from pymacruise.toads_dura import dura_initial_level_sets, run_dura_removal  # noqa: E402
from pymacruise.validation import (  # noqa: E402
    dice_coefficient,
    mesh_surface_distance_metrics,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--java-output-dir", type=Path, required=True)
    parser.add_argument("--base-stage", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    input_names = ("t1", "white_matter", "gray_matter", "white_matter_mask")
    images = {name: nib.load(str(args.input_dir / f"{name}.nii.gz")) for name in input_names}
    reference = images["t1"]
    if len(reference.shape) != 3 or tuple(reference.shape) != (172, 216, 164):
        raise ValueError("Dura Java comparison requires the 172x216x164 brain-resolution grid")
    if nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("Dura Java comparison inputs must be RAS-oriented")
    for name, image in images.items():
        if image.shape != reference.shape or not np.array_equal(image.affine, reference.affine):
            raise ValueError(f"{name} grid must exactly match the T1 grid")

    result = run_dura_removal(
        images["t1"], images["white_matter"], images["gray_matter"],
        images["white_matter_mask"], max_iterations=1, record_passes=True,
    )
    wm = np.asarray(result.wm_topology.dataobj, dtype=np.float32)
    gm_topology = np.asarray(result.gm_topology.dataobj, dtype=np.float32)
    gm = np.asarray(images["gray_matter"].dataobj, dtype=np.float32)
    initial = dura_initial_level_sets(wm, gm_topology)
    python_outputs = {
        "wm_tc": wm,
        "union": np.add(wm, gm, dtype=np.float32),
        "gm_tc": gm_topology,
        "inner_distance": initial.inner_signed_distance,
        "outer_distance": initial.outer_signed_distance,
        "level_set": np.asarray(result.evolution.level_set, dtype=np.float32),
        "stripped_t1": np.asarray(result.stripped_t1.dataobj, dtype=np.float32),
        "corrected_wm": np.asarray(result.corrected_white_matter.dataobj, dtype=np.float32),
        "corrected_gm": np.asarray(result.corrected_gray_matter.dataobj, dtype=np.float32),
        "brain_mask_level_set": np.asarray(result.brain_mask_level_set.dataobj, dtype=np.float32),
        "brain_mask": np.asarray(result.brain_mask.dataobj, dtype=np.uint8),
        "tgdm_pass_01": np.asarray(result.evolution.pass_level_sets[0], dtype=np.float32),
    }
    base_names = {
        "wm_tc": "wm_topology.nii.gz",
        "gm_tc": "gm_topology.nii.gz",
        "level_set": "brain_mask_level_set.nii.gz",
        "brain_mask_level_set": "brain_mask_level_set.nii.gz",
        "brain_mask": "brain_mask.nii.gz",
        "stripped_t1": "stripped_t1.nii.gz",
        "corrected_wm": "corrected_white_matter.nii.gz",
        "corrected_gm": "corrected_gray_matter.nii.gz",
    }
    outputs: dict[str, dict[str, object]] = {}
    for name, candidate in python_outputs.items():
        java = read_java_volume(
            args.java_output_dir / f"{name}.raw.gz", byte=name == "brain_mask"
        )
        error = np.abs(candidate.astype(np.float64) - java.astype(np.float64))
        item: dict[str, object] = {
            "exact": bool(np.array_equal(candidate, java)),
            "changed_voxels": int(np.count_nonzero(candidate != java)),
            "max_absolute_error": float(error.max(initial=0)),
            "mean_absolute_error": float(error.mean()),
        }
        if name in {"inner_distance", "outer_distance", "level_set", "brain_mask_level_set", "tgdm_pass_01"}:
            item["sign_changed_voxels"] = int(
                np.count_nonzero(np.signbit(candidate) != np.signbit(java))
            )
        if name == "brain_mask":
            item["dice_foreground"] = dice_coefficient(candidate > 0, java > 0)
        elif name in {"wm_tc", "gm_tc", "level_set", "tgdm_pass_01"}:
            item["dice_at_zero"] = dice_coefficient(candidate <= 0, java <= 0)
        if args.base_stage is not None and name in base_names:
            saved = np.asarray(nib.load(str(args.base_stage / base_names[name])).dataobj)
            item["python_recomputed_matches_workflow_stage"] = bool(
                np.array_equal(candidate, saved)
            )
        outputs[name] = item

    python_mask = python_outputs["brain_mask"] > 0
    java_mask = read_java_volume(args.java_output_dir / "brain_mask.raw.gz", byte=True) > 0

    def mask_mesh(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        vertices, faces, _, _ = marching_cubes(mask.astype(np.float32), level=0.5)
        return nib.affines.apply_affine(reference.affine, vertices), faces

    python_vertices, python_faces = mask_mesh(python_mask)
    java_vertices, java_faces = mask_mesh(java_mask)
    report = {
        "synthetic_only": True,
        "scope": "selected Java Dura wrapper on held-fixed brain-resolution RAS inputs",
        "shape": list(reference.shape),
        "orientation": "RAS",
        "spacing_mm": nib.affines.voxel_sizes(reference.affine).tolist(),
        "max_iterations": 1,
        "inputs_sha256": {
            name: hashlib.sha256((args.input_dir / f"{name}.nii.gz").read_bytes()).hexdigest()
            for name in input_names
        },
        "outputs": outputs,
        "mask_surface_metrics": mesh_surface_distance_metrics(
            python_vertices, python_faces, java_vertices, java_faces
        ),
    }
    if args.base_stage is not None:
        report["base_stage"] = str(args.base_stage)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "report": str(args.report),
        "brain_mask_dice": outputs["brain_mask"]["dice_foreground"],
        "mask_surface_metrics": report["mask_surface_metrics"],
        "outputs_compared": len(outputs),
    }, indent=2))


if __name__ == "__main__":
    main()
