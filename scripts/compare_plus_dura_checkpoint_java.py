#!/usr/bin/env python3
"""Compare a saved PLUS Dura checkpoint with selected Java source loops.

This replays the selected ``LesionToadSegmentation.generateDuraRemovalOutputs``
method from the exact registered-loop memberships and segmentation saved by a
PLUS workflow. It applies the wrapper's subsequent float32 scale by 255 before
comparing with the saved Python Dura products. This isolates that adapter; it
does not compare upstream TOADS/LesionTOADS evolution or native-grid transfer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np


REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "tests/java_reference/LesionToadsDuraCheckpointReferenceHarness.java"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CHANNELS = ("csf", "gm", "wm_fill", "wm_mask")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checked_input(record: dict, label: str) -> Path:
    path = Path(record["path"])
    if not path.is_file() or sha256(path) != record.get("sha256"):
        raise ValueError(f"{label} missing or failed SHA-256 check: {path}")
    return path


def compare(workflow_root: Path, output_dir: Path, java_xmx: str) -> dict:
    workflow_root = workflow_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = workflow_root / "dura-outputs/lesion-dura-output-checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text())
    if checkpoint.get("stage") != "plus-lesion-dura-output-checkpoint":
        raise ValueError("workflow has no completed PLUS Dura checkpoint")
    inputs = checkpoint["inputs"]
    memberships_path = checked_input(inputs["registered_memberships"], "registered memberships")
    segmentation_path = checked_input(inputs["registered_segmentation"], "registered segmentation")
    names = checkpoint["class_names"]
    if len(names) != 13:
        raise ValueError(f"expected 13 TOADS class names; found {len(names)}")
    membership_image = nib.load(str(memberships_path))
    segmentation_image = nib.load(str(segmentation_path))
    shape = segmentation_image.shape
    if membership_image.shape != (*shape, 17):
        raise ValueError("registered memberships must contain 13 classes, lesion, and black-hole channels")
    if nib.aff2axcodes(segmentation_image.affine) != ("R", "A", "S"):
        raise ValueError("registered segmentation must be RAS-oriented")
    if not np.allclose(membership_image.affine, segmentation_image.affine, rtol=0, atol=1e-5):
        raise ValueError("registered memberships and segmentation are on different grids")

    memberships = np.asarray(membership_image.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_image.dataobj, dtype=np.int32)
    input_path = output_dir / "java-input.bin"
    with input_path.open("wb") as stream:
        stream.write(struct.pack(">IIII", *shape, len(names)))
        stream.write(memberships[..., :13].astype(">f4", copy=False).tobytes(order="C"))
        stream.write(segmentation.astype(">i4", copy=False).tobytes(order="C"))
        stream.write(memberships[..., 13].astype(">f4", copy=False).tobytes(order="C"))
        stream.write(memberships[..., 14].astype(">f4", copy=False).tobytes(order="C"))
    del memberships, segmentation

    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("JDK javac and java are required")
    classes = output_dir / "classes"
    classes.mkdir(exist_ok=True)
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True)
    java_output = output_dir / "java-output.bin"
    subprocess.run(
        [java, f"-Xms256m", f"-Xmx{java_xmx}", "-cp", str(classes),
         HARNESS.stem, str(input_path), str(java_output), *names],
        check=True,
    )
    voxel_count = int(np.prod(shape))
    source_outputs = np.memmap(java_output, dtype=">f4", mode="r", shape=(4, voxel_count))
    records = {}
    for index, name in enumerate(CHANNELS):
        path = checked_input(checkpoint["outputs"][name], f"Python {name}")
        image = nib.load(str(path))
        if image.shape != shape or not np.allclose(image.affine, segmentation_image.affine, rtol=0, atol=1e-5):
            raise ValueError(f"Python {name} output does not share the RAS checkpoint grid")
        java_scaled = np.asarray(source_outputs[index], dtype=np.float32) * np.float32(255.0)
        python = np.asarray(image.dataobj, dtype=np.float32).ravel(order="C")
        delta = np.abs(java_scaled.astype(np.float64) - python.astype(np.float64))
        records[name] = {
            "exact_float32_after_source_scale": bool(np.array_equal(java_scaled, python)),
            "changed_voxels": int(np.count_nonzero(java_scaled != python)),
            "voxel_count": voxel_count,
            "max_absolute_error": float(delta.max(initial=0.0)),
            "mean_absolute_error": float(delta.mean()),
        }
    return {
        "case": "saved corrected atlas-aligned synthetic PLUS run",
        "synthetic_only": True,
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "source_method": "LesionToadSegmentation.generateDuraRemovalOutputs followed by scaleFloatValue(255)",
        "selected_java_source_sha256": SOURCE_SHA256,
        "checkpoint_run_signature": checkpoint["run_signature"],
        "inputs": {key: value["sha256"] for key, value in inputs.items() if key in {"registered_memberships", "registered_segmentation"}},
        "outputs": records,
        "all_exact": all(record["exact_float32_after_source_scale"] for record in records.values()),
        "scope": "Dura output adapter only; upstream lesion memberships and native-grid transfer are not validated by this comparison",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--java-xmx", default="6g", help="maximum Java heap (default: 6g)")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.workflow_root, args.output_dir, args.java_xmx)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_exact"]:
        raise SystemExit("selected Java and saved Python Dura outputs differ")


if __name__ == "__main__":
    main()
