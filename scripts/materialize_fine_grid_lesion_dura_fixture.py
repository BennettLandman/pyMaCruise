#!/usr/bin/env python3
"""Preserve the exact P41z iteration-6 state used for a Java Dura comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import nifti_ras_on_reference, save_nifti_ras


REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / "build/plus-vandy-atlas-aligned-fine"
DEST = REPO / "tests/data/synthetic/plus_fine_grid_dura_iteration_6"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def checked(record: dict, label: str) -> Path:
    path = Path(record["path"]).resolve()
    if not path.is_file() or path.is_symlink() or digest(path) != record.get("sha256"):
        raise ValueError(f"{label} missing or fails checkpoint hash: {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEST)
    args = parser.parse_args()
    root = args.workflow_root.resolve()
    dest = args.output_dir.resolve()
    loop_manifest_path = root / "loop-chunks/iteration-0006/lesion-registered-loop.json"
    dura_manifest_path = root / "surface-checkpoint-0006/lesion-dura/lesion-dura-output-checkpoint.json"
    loop = json.loads(loop_manifest_path.read_text())
    dura = json.loads(dura_manifest_path.read_text())
    if loop.get("status") != "complete" or dura.get("status") != "complete":
        raise ValueError("both source checkpoints must be complete")
    summary_path = checked(loop["outputs"]["summary"], "registered loop summary")
    summary = json.loads(summary_path.read_text())
    iteration = int(summary["completed_iteration"])
    if iteration != 6:
        raise ValueError(f"expected checkpoint completed iteration 6, got {iteration}")
    loop_inputs = loop["inputs"]
    recrop_manifest_path = checked(loop_inputs["registered_recrop_manifest"], "registered recrop manifest")
    recrop_manifest = json.loads(recrop_manifest_path.read_text())
    recrop_summary_path = checked(recrop_manifest["outputs"]["summary"], "registered recrop summary")
    class_names = json.loads(recrop_summary_path.read_text())["class_names"]
    if len(class_names) != 13:
        raise ValueError("the saved loop must have exactly 13 atlas memberships")
    memberships_path = checked(loop["outputs"]["memberships"], "registered memberships")
    segmentation_path = checked(loop["outputs"]["segmentation"], "registered segmentation")
    membership_image = nib.load(str(memberships_path))
    segmentation_image = nib.load(str(segmentation_path))
    if membership_image.shape != (*segmentation_image.shape, 17):
        raise ValueError("registered membership checkpoint must have 17 channels")
    if nib.aff2axcodes(membership_image.affine) != ("R", "A", "S"):
        raise ValueError("checkpoint memberships are not in RAS")
    if not np.allclose(membership_image.affine, segmentation_image.affine, rtol=0, atol=1e-5):
        raise ValueError("checkpoint memberships and segmentation grid differ")
    packed = np.asarray(membership_image.dataobj, dtype=np.float32)
    dest.mkdir(parents=True, exist_ok=True)
    # Keep only the 13 class probabilities and the two exact auxiliary maps
    # consumed by generateDuraRemovalOutputs; other checkpoint channels are
    # not silently substituted into this source comparison.
    saved = {}
    for name, values in (
        ("memberships_15ch", packed[..., :15]),
        ("segmentation", np.asarray(segmentation_image.dataobj, dtype=np.int32)),
    ):
        path = dest / f"{name}.nii.gz"
        save_nifti_ras(nifti_ras_on_reference(segmentation_image, values), path)
        saved[name] = {"path": path.name, "sha256": digest(path)}
    for name in ("csf", "gm", "wm_fill", "wm_mask"):
        source_path = checked(dura["outputs"][name], f"saved Python {name}")
        image = nib.load(str(source_path))
        if image.shape != segmentation_image.shape or not np.allclose(image.affine, segmentation_image.affine, rtol=0, atol=1e-5):
            raise ValueError(f"saved Python output {name} has a mismatched grid")
        target = dest / f"lesion_toads_{name}.nii.gz"
        # Re-save the exact float32 voxel array through the RAS helper; the
        # source checkpoint is already RAS and this preserves its voxel bits.
        values = np.asarray(image.dataobj, dtype=np.float32)
        save_nifti_ras(nifti_ras_on_reference(segmentation_image, values), target, dtype=np.float32)
        saved[f"python_{name}"] = {"path": target.name, "sha256": digest(target)}
    provenance = {
        "case": "P41z fine-grid synthetic LesionTOADS Dura outputs",
        "selected_java_source": "VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java:generateDuraRemovalOutputs",
        "selected_java_source_sha256": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
        "synthetic_only": True,
        "completed_iteration": iteration,
        "energy_ratio": summary.get("energy_ratio"),
        "spatial_convention": "RAS",
        "shape_xyz": list(segmentation_image.shape),
        "class_names": class_names,
        "membership_layout": {
            "classes": [0, 12],
            "lesion_membership": 13,
            "black_hole_membership": 14,
            "only_channels_saved": 15,
        },
        "source_checkpoint": {
            "loop_manifest_sha256": digest(loop_manifest_path),
            "dura_manifest_sha256": digest(dura_manifest_path),
            "registered_memberships_sha256": digest(memberships_path),
            "registered_segmentation_sha256": digest(segmentation_path),
            "registered_recrop_manifest_sha256": digest(recrop_manifest_path),
            "registered_recrop_summary_sha256": digest(recrop_summary_path),
        },
        "files": saved,
    }
    manifest_path = dest / "provenance.json"
    manifest_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"fixture_dir": str(dest), "provenance": str(manifest_path),
                      "image_count": len(saved), "size_bytes": sum((dest / item["path"]).stat().st_size for item in saved.values())}, indent=2))


if __name__ == "__main__":
    main()
