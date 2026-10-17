#!/usr/bin/env python3
"""Rebuild exact held-fixed inputs for a saved plus-workflow TGDM stage."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import struct

import nibabel as nib
import numpy as np

from pymacruise.plus_pipeline import (
    PlusPipelineOptions,
    run_plus_ace_stage,
    run_plus_boundary_stage,
    run_plus_gvf_stage,
    run_plus_tissue_update_stage,
)
from pymacruise.plus_topology_composition import run_plus_topology_and_tgdm
from pymacruise.boundary import multi_atlas_csf_boundary_image
from pymacruise.ace import _tgdm_signed_distance_fmm
from pymacruise.tgdm import _SMOOTHING_PROFILE_SETTINGS, tgdm_evolve

INPUT_ORDER = (
    "segmentation", "toads_wm_mask", "multi_atlas_gm", "multi_atlas_wm",
    "multi_atlas_mask", "toads_gm", "toads_wm_fill", "toads_wm_mask",
)
OUTPUT_NAMES = ("initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True,
                        help="saved plus-surface workflow directory containing run.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--context-manifest",
        type=Path,
        help="optional composition manifest when the workflow is checkpoint-composed",
    )
    parser.add_argument(
        "--use-saved-tgdm-outputs",
        action="store_true",
        help="reuse saved TGDM products and rebuild only held-fixed upstream arrays",
    )
    parser.add_argument(
        "--allow-python-stage-differences",
        action="store_true",
        help="record, rather than abort on, replay differences from saved Python TGDM products",
    )
    args = parser.parse_args()
    root = args.workflow_root.resolve()
    manifest_path = root / "run.json"
    manifest = json.loads(manifest_path.read_text())
    crop_outputs = manifest["stages"]["transfer_crop"]["outputs"]
    inputs = {
        name: nib.load(crop_outputs[name]["path"])
        for name in set(INPUT_ORDER) | {"brain_mask"}
    }
    settings = dict(manifest["settings"]["plus_pipeline"])
    settings["boundary_overlook_labels"] = tuple(settings["boundary_overlook_labels"])
    options = PlusPipelineOptions(**settings)
    saved_stage = root / "stages/topology_ace_gvf_tgdm" / manifest["stages"]["topology_ace_gvf_tgdm"]["attempt"]
    saved_stage_outputs = manifest["stages"]["topology_ace_gvf_tgdm"]["outputs"]
    python_outputs: dict[str, np.ndarray] = {}
    saved_stage_differences: dict[str, int] = {}
    if args.use_saved_tgdm_outputs:
        # Rebuild only TGDM's held-fixed input fields. The terminal stage has
        # already persisted the five output arrays, and repeating nested TGDM
        # on this grid is a multi-minute operation.
        boundary = run_plus_boundary_stage(inputs["segmentation"], inputs["toads_wm_fill"], options)
        tissue = run_plus_tissue_update_stage(*(inputs[name] for name in INPUT_ORDER[2:]))
        ace = run_plus_ace_stage(tissue.gm, tissue.wm, boundary.cleaned_boundary, options)
        gvf = run_plus_gvf_stage(ace.enhanced_gray_matter, options)
        for name in OUTPUT_NAMES:
            record = saved_stage_outputs[name]
            path = Path(record["path"])
            if not path.is_file() or sha256(path) != record["sha256"]:
                raise SystemExit(f"saved TGDM output hash check failed: {name}")
            python_outputs[name] = np.asarray(nib.load(str(path)).dataobj, dtype=np.float32)
        tgdm_initial_membership = python_outputs["initial_level_set"]
        tissue_wm, tissue_gm = tissue.wm, tissue.gm
        ace_gm, ace_skeleton, gvf_field = ace.enhanced_gray_matter, ace.thinned_skeleton, gvf.vector_field
        ace_full_skeleton = ace.skeleton
        java_input_mode = "precomputed_initial_level_set"
    else:
        result = run_plus_topology_and_tgdm(
            *(inputs[name] for name in INPUT_ORDER),
            options=options,
            brain_mask=inputs["brain_mask"],
        )
        for name in OUTPUT_NAMES:
            actual = np.asarray(getattr(result.tgdm, name).dataobj)
            saved = np.asarray(nib.load(saved_stage / f"{name}.nii.gz").dataobj)
            if not np.array_equal(actual, saved):
                changed = int(np.count_nonzero(actual != saved))
                saved_stage_differences[name] = changed
                if not args.allow_python_stage_differences:
                    raise SystemExit(f"Python replay differs from saved workflow stage {name}: {changed} voxels")
            python_outputs[name] = actual.astype(np.float32, copy=False)
        tgdm_initial_membership = result.tgdm_initial_membership
        tissue_wm, tissue_gm = result.tissue.wm, result.tissue.gm
        ace_gm, ace_skeleton, gvf_field = (
            result.ace.enhanced_gray_matter,
            result.ace.thinned_skeleton,
            result.gvf.vector_field,
        )
        ace_full_skeleton = result.ace.skeleton
        boundary = result.boundary
        smoothed_initial = result.tgdm.execution.smoothed_initial.level_set
        java_input_mode = "initial_membership"
    if args.use_saved_tgdm_outputs:
        # The workflow manifest intentionally persists only delivered TGDM
        # products. Rebuild the intermediate smoothing surface separately so
        # source comparisons can localize where a later surface mismatch starts.
        smoothing_input = _tgdm_signed_distance_fmm(tgdm_initial_membership, max_distance=8.0)
        smoothed_initial = tgdm_evolve(
            smoothing_input,
            profile="central",
            curvature_force=-1.0,
            pressure_force=0.0,
            outer_iterations=6,
            settings=_SMOOTHING_PROFILE_SETTINGS,
        ).level_set
    held = (
        np.asarray(tgdm_initial_membership.dataobj if hasattr(tgdm_initial_membership, "dataobj") else tgdm_initial_membership, dtype=np.float32),
        np.asarray(tissue_wm.dataobj, dtype=np.float32),
        np.asarray(tissue_gm.dataobj, dtype=np.float32),
        np.asarray(ace_gm.dataobj, dtype=np.float32),
        np.asarray(ace_skeleton.dataobj, dtype=np.float32),
        np.asarray(gvf_field.dataobj, dtype=np.float32),
    )
    shape = held[0].shape
    if any(array.shape != shape for array in held[:5]) or held[5].shape != (*shape, 3):
        raise SystemExit("reconstructed plus TGDM inputs do not share the expected scalar/vector grids")
    if any(not np.isfinite(array).all() for array in held):
        raise SystemExit("reconstructed plus TGDM inputs contain non-finite values")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_path = args.output_dir / "arrays.bin"
    with input_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for array in held:
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="C"))
    np.save(args.output_dir / "python-held-initial-level-set.npy", held[0])
    binary_brain_mask = np.asarray(inputs["brain_mask"].dataobj, dtype=np.float32)
    if binary_brain_mask.shape != shape or not np.isfinite(binary_brain_mask).all():
        raise SystemExit("cropped brain mask must be finite and match the TGDM input grid")
    brain_mask = np.where(binary_brain_mask > 0, -1.0, 1.0).astype(np.float32)
    np.save(args.output_dir / "python-brain-mask-binary.npy", binary_brain_mask)
    brain_mask_path = args.output_dir / "brain-mask.bin"
    with brain_mask_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        stream.write(brain_mask.astype(">f4", copy=False).tobytes(order="C"))
    np.save(args.output_dir / "python-brain-mask-bound.npy", brain_mask)
    ace_inputs=(np.asarray(tissue_gm.dataobj,dtype=np.float32),np.asarray(tissue_wm.dataobj,dtype=np.float32),np.asarray(boundary.cleaned_boundary.dataobj,dtype=np.float32))
    ace_input_path=args.output_dir/"ace-input.bin"
    with ace_input_path.open("wb") as stream:
        stream.write(struct.pack(">III",*shape))
        for array in ace_inputs: stream.write(np.asarray(array,dtype=">f4").tobytes(order="C"))
    segmentation = np.asarray(inputs["segmentation"].dataobj, dtype=np.int32)
    raw_boundary = np.asarray(boundary.raw_boundary.dataobj, dtype=np.int32)
    boundary_input_path = args.output_dir / "boundary-labels.bin"
    with boundary_input_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        stream.write(segmentation.astype(">i4", copy=False).tobytes(order="C"))
    np.save(args.output_dir / "python-raw-boundary.npy", raw_boundary)
    folded_fixture_path = Path(__file__).resolve().parents[1] / "tests/data/synthetic/plus_boundary_folded_cortex/segmentation.nii.gz"
    folded_fixture = nib.load(str(folded_fixture_path))
    folded_fixture_compatible = (
        folded_fixture.shape == shape
        and np.allclose(folded_fixture.affine, nib.load(str(saved_stage_outputs["inner"]["path"])).affine, rtol=0, atol=1e-5)
    )
    folded_boundary_input_sha256 = None
    folded_raw_boundary_sha256 = None
    if folded_fixture_compatible:
        folded_boundary = np.asarray(multi_atlas_csf_boundary_image(folded_fixture).dataobj, dtype=np.int32)
        folded_segmentation = np.asarray(folded_fixture.dataobj, dtype=np.int32)
        folded_input_path = args.output_dir / "boundary-folded-labels.bin"
        with folded_input_path.open("wb") as stream:
            stream.write(struct.pack(">III", *shape))
            stream.write(folded_segmentation.astype(">i4", copy=False).tobytes(order="C"))
        np.save(args.output_dir / "python-folded-raw-boundary.npy", folded_boundary)
        folded_boundary_input_sha256 = sha256(folded_input_path)
        folded_raw_boundary_sha256 = hashlib.sha256(
            folded_boundary.astype(">i4", copy=False).tobytes(order="C")
        ).hexdigest()
    gvf_input_path=args.output_dir/"gvf-input.bin"
    with gvf_input_path.open("wb") as stream:
        stream.write(struct.pack(">III",*shape))
        stream.write(np.asarray(ace_gm.dataobj,dtype=">f4").tobytes(order="C"))
    for name,array in (("ace-enhanced",ace_gm),("ace-skeleton",ace_full_skeleton),("ace-thinned-skeleton",ace_skeleton)):
        np.save(args.output_dir/f"python-{name}.npy",np.asarray(array.dataobj))
    for name,array in (("tgdm-wm",tissue_wm),("tgdm-gm",tissue_gm)):
        np.save(args.output_dir/f"python-{name}.npy",np.asarray(array.dataobj,dtype=np.float32))
    np.save(args.output_dir/"python-gvf-vector.npy",np.asarray(gvf_field.dataobj,dtype=np.float32))
    for name in OUTPUT_NAMES:
        np.save(args.output_dir / f"python-{name}.npy", python_outputs[name])
    np.save(args.output_dir / "python-smoothed_initial.npy", smoothed_initial)
    full_manifest_path = (
        args.context_manifest.resolve()
        if args.context_manifest is not None
        else root.parent / "plus-full-workflow.json"
    )
    full_manifest = json.loads(full_manifest_path.read_text())
    primary_input_paths = [
        Path(record["path"]).resolve()
        for name, record in full_manifest.get("inputs", {}).items()
        if name in {"t1", "flair", "brain_mask", "maseg"} and record.get("path")
    ]
    synthetic_root = Path(__file__).resolve().parents[1] / "tests/data/synthetic"
    synthetic_only = bool(primary_input_paths) and all(
        path.is_relative_to(synthetic_root) for path in primary_input_paths
    )
    evidence_category = (
        "synthetic held-fixed source comparison"
        if synthetic_only
        else "real-data held-fixed source comparison"
    )
    input_records = {
        name: {"path": crop_outputs[name]["path"], "sha256": sha256(Path(crop_outputs[name]["path"]))}
        for name in sorted(set(INPUT_ORDER))
    }
    input_records["brain_mask"] = {
        "path": crop_outputs["brain_mask"]["path"],
        "sha256": sha256(Path(crop_outputs["brain_mask"]["path"])),
    }
    report = {
        "synthetic_only": synthetic_only,
        "evidence_category": evidence_category,
        "workflow_manifest": str(manifest_path),
        "workflow_manifest_sha256": sha256(manifest_path),
        "workflow_run_signature": manifest["run_signature"],
        "full_workflow_manifest": str(full_manifest_path),
        "full_workflow_manifest_sha256": sha256(full_manifest_path),
        "full_workflow_run_signature": full_manifest["run_signature"],
        "implementation_signature": full_manifest.get("implementation_signature"),
        "context_manifest_kind": (
            "checkpoint_composition" if args.context_manifest is not None else "plus_full_workflow"
        ),
        "tgdm_input_mode": java_input_mode,
        "shape_xyz": list(shape),
        "affine_ras": np.asarray(nib.load(str(saved_stage_outputs["inner"]["path"])).affine).tolist(),
        "settings": options.as_dict(),
        "crop_inputs": input_records,
        "replayed_saved_outputs_exactly": (
            [] if args.use_saved_tgdm_outputs else [name for name in OUTPUT_NAMES if name not in saved_stage_differences]
        ),
        "saved_stage_difference_voxels": saved_stage_differences,
        "ace_gvf_parameters": {"csf_merge_mode":3,"gvf_lambda":0.2,"gvf_iterations":200,"gvf_normalize":False},
        "held_fixed_input_binary_sha256": sha256(input_path),
        "brain_mask_bound_binary_sha256": sha256(brain_mask_path),
        "brain_mask_semantics": "Python adapter's signed containment bound: inside=-1, outside=+1; source JIST plugin input is explicitly named Brain Mask Level Set",
        "ace_input_binary_sha256": sha256(ace_input_path),
        "boundary_input_binary_sha256": sha256(boundary_input_path),
        "raw_boundary_sha256": hashlib.sha256(raw_boundary.astype(">i4", copy=False).tobytes(order="C")).hexdigest(),
        "folded_boundary_fixture_compatible": folded_fixture_compatible,
        "folded_boundary_input_sha256": folded_boundary_input_sha256,
        "folded_raw_boundary_sha256": folded_raw_boundary_sha256,
        "gvf_input_binary_sha256": sha256(gvf_input_path),
    }
    (args.output_dir / "input-provenance.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
