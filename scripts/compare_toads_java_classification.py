#!/usr/bin/env python3
"""Compare a source Java TOADS byte map with a Python NIfTI classification.

The Java harness exports ModelImage's x-fastest linear buffer. For the default
RAS NIfTI array convention, reshape its bytes with Fortran order.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def compare(java_raw: Path, python_nifti: Path, *, order: str = "F") -> dict[str, object]:
    image = nib.load(str(python_nifti))
    if len(image.shape) != 3:
        raise ValueError("classification NIfTI must be 3-D")
    python = np.asarray(image.dataobj).astype(np.uint8, copy=False)
    java_bytes = np.fromfile(java_raw, dtype=np.uint8)
    if java_bytes.size != python.size:
        raise ValueError(
            f"Java byte count {java_bytes.size} does not match NIfTI voxel count {python.size}"
        )
    java = java_bytes.reshape(python.shape, order=order)
    labels = sorted(set(np.unique(python).tolist()) | set(np.unique(java).tolist()))
    by_label: dict[str, object] = {}
    for label in labels:
        p = python == label
        j = java == label
        denominator = int(p.sum() + j.sum())
        by_label[str(label)] = {
            "python_voxels": int(p.sum()),
            "java_voxels": int(j.sum()),
            "dice": float(2 * np.count_nonzero(p & j) / denominator) if denominator else 1.0,
        }
    return {
        "python_nifti": str(python_nifti),
        "java_raw": str(java_raw),
        "java_reshape_order": order,
        "shape": list(python.shape),
        "voxel_count": int(python.size),
        "changed_voxels": int(np.count_nonzero(python != java)),
        "labels": by_label,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("java_raw", type=Path)
    parser.add_argument("python_nifti", type=Path)
    parser.add_argument("--java-reshape-order", choices=("C", "F"), default="F")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = compare(args.java_raw, args.python_nifti, order=args.java_reshape_order)
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
