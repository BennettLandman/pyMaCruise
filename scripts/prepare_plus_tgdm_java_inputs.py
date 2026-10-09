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

from pymacruise.plus_pipeline import PlusPipelineOptions
from pymacruise.plus_topology_composition import run_plus_topology_and_tgdm
from pymacruise.boundary import multi_atlas_csf_boundary_image

INPUT_ORDER = (
    "segmentation", "toads_wm_fill", "multi_atlas_gm", "multi_atlas_wm",
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
    args = parser.parse_args()
    root = args.workflow_root.resolve()
    manifest_path = root / "run.json"
    manifest = json.loads(manifest_path.read_text())
    crop_outputs = manifest["stages"]["transfer_crop"]["outputs"]
    inputs = {name: nib.load(crop_outputs[name]["path"]) for name in set(INPUT_ORDER)}
    settings = dict(manifest["settings"]["plus_pipeline"])
    settings["boundary_overlook_labels"] = tuple(settings["boundary_overlook_labels"])
    options = PlusPipelineOptions(**settings)
    result = run_plus_topology_and_tgdm(*(inputs[name] for name in INPUT_ORDER), options=options)

    saved_stage = root / "stages/topology_ace_gvf_tgdm" / manifest["stages"]["topology_ace_gvf_tgdm"]["attempt"]
    for name in OUTPUT_NAMES:
        actual = np.asarray(getattr(result.tgdm, name).dataobj)
        saved = np.asarray(nib.load(saved_stage / f"{name}.nii.gz").dataobj)
        if not np.array_equal(actual, saved):
            changed = int(np.count_nonzero(actual != saved))
            raise SystemExit(f"Python replay differs from saved workflow stage {name}: {changed} voxels")

    held = (
        np.asarray(result.tgdm_initial_membership.dataobj, dtype=np.float32),
        np.asarray(result.tissue.wm.dataobj, dtype=np.float32),
        np.asarray(result.tissue.gm.dataobj, dtype=np.float32),
        np.asarray(result.ace.enhanced_gray_matter.dataobj, dtype=np.float32),
        np.asarray(result.ace.thinned_skeleton.dataobj, dtype=np.float32),
        np.asarray(result.gvf.vector_field.dataobj, dtype=np.float32),
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
    ace_inputs=(np.asarray(result.tissue.gm.dataobj,dtype=np.float32),np.asarray(result.tissue.wm.dataobj,dtype=np.float32),np.asarray(result.boundary.cleaned_boundary.dataobj,dtype=np.float32))
    ace_input_path=args.output_dir/"ace-input.bin"
    with ace_input_path.open("wb") as stream:
        stream.write(struct.pack(">III",*shape))
        for array in ace_inputs: stream.write(np.asarray(array,dtype=">f4").tobytes(order="C"))
    segmentation = np.asarray(inputs["segmentation"].dataobj, dtype=np.int32)
    raw_boundary = np.asarray(result.boundary.raw_boundary.dataobj, dtype=np.int32)
    boundary_input_path = args.output_dir / "boundary-labels.bin"
    with boundary_input_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        stream.write(segmentation.astype(">i4", copy=False).tobytes(order="C"))
    np.save(args.output_dir / "python-raw-boundary.npy", raw_boundary)
    folded_fixture_path = Path(__file__).resolve().parents[1] / "tests/data/synthetic/plus_boundary_folded_cortex/segmentation.nii.gz"
    folded_fixture = nib.load(str(folded_fixture_path))
    if folded_fixture.shape != shape or not np.allclose(folded_fixture.affine, result.tgdm.inner.affine, rtol=0, atol=1e-5):
        raise SystemExit("folded boundary fixture does not match the held-fixed plus crop grid")
    folded_boundary = np.asarray(multi_atlas_csf_boundary_image(folded_fixture).dataobj, dtype=np.int32)
    folded_segmentation = np.asarray(folded_fixture.dataobj, dtype=np.int32)
    folded_input_path = args.output_dir / "boundary-folded-labels.bin"
    with folded_input_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        stream.write(folded_segmentation.astype(">i4", copy=False).tobytes(order="C"))
    np.save(args.output_dir / "python-folded-raw-boundary.npy", folded_boundary)
    gvf_input_path=args.output_dir/"gvf-input.bin"
    with gvf_input_path.open("wb") as stream:
        stream.write(struct.pack(">III",*shape))
        stream.write(np.asarray(result.ace.enhanced_gray_matter.dataobj,dtype=">f4").tobytes(order="C"))
    for name,array in (("ace-enhanced",result.ace.enhanced_gray_matter),("ace-skeleton",result.ace.skeleton),("ace-thinned-skeleton",result.ace.thinned_skeleton)):
        np.save(args.output_dir/f"python-{name}.npy",np.asarray(array.dataobj))
    np.save(args.output_dir/"python-gvf-vector.npy",np.asarray(result.gvf.vector_field.dataobj,dtype=np.float32))
    for name in OUTPUT_NAMES:
        np.save(args.output_dir / f"python-{name}.npy", np.asarray(getattr(result.tgdm, name).dataobj))
    full_manifest_path = root.parent / "plus-full-workflow.json"
    full_manifest = json.loads(full_manifest_path.read_text())
    input_records = {
        name: {"path": crop_outputs[name]["path"], "sha256": sha256(Path(crop_outputs[name]["path"]))}
        for name in sorted(set(INPUT_ORDER))
    }
    report = {
        "synthetic_only": True,
        "workflow_manifest": str(manifest_path),
        "workflow_manifest_sha256": sha256(manifest_path),
        "workflow_run_signature": manifest["run_signature"],
        "full_workflow_manifest": str(full_manifest_path),
        "full_workflow_manifest_sha256": sha256(full_manifest_path),
        "full_workflow_run_signature": full_manifest["run_signature"],
        "implementation_signature": full_manifest["implementation_signature"],
        "shape_xyz": list(shape),
        "affine_ras": np.asarray(result.tgdm.inner.affine).tolist(),
        "settings": options.as_dict(),
        "crop_inputs": input_records,
        "replayed_saved_outputs_exactly": list(OUTPUT_NAMES),
        "ace_gvf_parameters": {"csf_merge_mode":3,"gvf_lambda":0.2,"gvf_iterations":200,"gvf_normalize":False},
        "held_fixed_input_binary_sha256": sha256(input_path),
        "ace_input_binary_sha256": sha256(ace_input_path),
        "boundary_input_binary_sha256": sha256(boundary_input_path),
        "raw_boundary_sha256": hashlib.sha256(raw_boundary.astype(">i4", copy=False).tobytes(order="C")).hexdigest(),
        "folded_boundary_input_sha256": sha256(folded_input_path),
        "folded_raw_boundary_sha256": hashlib.sha256(folded_boundary.astype(">i4", copy=False).tobytes(order="C")).hexdigest(),
        "gvf_input_binary_sha256": sha256(gvf_input_path),
    }
    (args.output_dir / "input-provenance.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
