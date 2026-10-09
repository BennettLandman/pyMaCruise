#!/usr/bin/env python3
"""Compare Java TOADS DuraRemoval_input products with saved Python NIfTIs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


FLOAT_PRODUCTS = ("csf", "gm", "wmfill", "wmmask")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(java_dir: Path, python_dir: Path, java_scale: float) -> dict[str, object]:
    segmentation_paths = sorted(python_dir.glob("*_Toads_seg.nii.gz"))
    if len(segmentation_paths) != 1:
        raise ValueError("Python output directory must contain exactly one *_Toads_seg.nii.gz")
    python_seg_path = segmentation_paths[0]
    prefix = python_seg_path.name.removesuffix("_Toads_seg.nii.gz")
    java_seg_path = java_dir / "classification.raw"
    python_seg_image = nib.load(str(python_seg_path))
    if nib.aff2axcodes(python_seg_image.affine) != ("R", "A", "S"):
        raise ValueError("Python classification must be RAS")
    shape = python_seg_image.shape
    python_seg = np.asarray(python_seg_image.dataobj, dtype=np.uint8)
    java_seg_raw = np.fromfile(java_seg_path, dtype=np.uint8)
    if java_seg_raw.size != python_seg.size:
        raise ValueError("Java and Python classifications have different voxel counts")
    java_seg = java_seg_raw.reshape(shape, order="F")
    outputs: dict[str, object] = {
        "classification": {
            "exact": bool(np.array_equal(java_seg, python_seg)),
            "changed_voxels": int(np.count_nonzero(java_seg != python_seg)),
            "voxel_count": int(python_seg.size),
            "java_sha256": sha256(java_seg_path),
            "python_sha256": sha256(python_seg_path),
        }
    }
    for index, name in enumerate(FLOAT_PRODUCTS, start=1):
        java_path = java_dir / f"{name}.raw"
        python_path = python_dir / f"{prefix}_{name}.nii.gz"
        python_image = nib.load(str(python_path))
        if python_image.shape != shape or nib.aff2axcodes(python_image.affine) != ("R", "A", "S"):
            raise ValueError(f"Python {name} image must use the classification RAS grid")
        if not np.allclose(python_image.affine, python_seg_image.affine, rtol=0, atol=1e-5):
            raise ValueError(f"Python {name} image affine differs from classification")
        java = np.fromfile(java_path, dtype="<f4")
        if java.size != python_seg.size:
            raise ValueError(f"Java {name} field has the wrong voxel count")
        java = java.reshape(shape, order="F").astype(np.float32, copy=False)
        python = np.asarray(python_image.dataobj, dtype=np.float32)
        scaled_java = np.multiply(java, np.float32(java_scale), dtype=np.float32)
        delta = np.abs(scaled_java.astype(np.float64) - python.astype(np.float64))
        outputs[name] = {
            "java_scale": java_scale,
            "exact_after_scale": bool(np.array_equal(scaled_java, python)),
            "changed_voxels_after_scale": int(np.count_nonzero(scaled_java != python)),
            "max_absolute_error": float(delta.max(initial=0.0)),
            "mean_absolute_error": float(delta.mean()),
            "java_sha256": sha256(java_path),
            "python_sha256": sha256(python_path),
        }
    return {
        "synthetic_only": True,
        "orientation": "RAS",
        "shape": list(shape),
        "java_output_mode": "DuraRemoval_input",
        "java_linear_order": "x-fastest / Fortran reshape",
        "java_unconditional_wrapper_scale": java_scale,
        "outputs": outputs,
        "all_outputs_exact": all(
            bool(value.get("exact", value.get("exact_after_scale", False)))
            for value in outputs.values()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("java_output_dir", type=Path)
    parser.add_argument("python_output_dir", type=Path)
    parser.add_argument("--java-scale", type=float, default=255.0)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = compare(args.java_output_dir, args.python_output_dir, args.java_scale)
    serialized = json.dumps(result, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n", encoding="utf-8")
    if not result["all_outputs_exact"]:
        raise SystemExit("Java and Python TOADS DuraRemoval_input outputs differ")


if __name__ == "__main__":
    main()
