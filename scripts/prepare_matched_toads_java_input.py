#!/usr/bin/env python3
"""Serialize a held-fixed masked RAS T1 for the selected Java TOADS harness.

Java's test harness accepts one intensity image and has no brain-mask argument.
This writes the same mask-applied intensity array that Python TOADS prepares,
leaving robust intensity normalization to each implementation. The array is
little-endian float32 in MIPAV x-fastest order.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t1-n4", required=True, type=Path)
    parser.add_argument("--brain-mask", required=True, type=Path)
    parser.add_argument("--raw-output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    t1_path, mask_path = args.t1_n4.resolve(), args.brain_mask.resolve()
    t1, mask = nib.load(str(t1_path)), nib.load(str(mask_path))
    if len(t1.shape) != 3 or len(mask.shape) != 3:
        raise ValueError("held-fixed T1 and mask must be three-dimensional")
    if nib.aff2axcodes(t1.affine) != ("R", "A", "S") or nib.aff2axcodes(mask.affine) != ("R", "A", "S"):
        raise ValueError("held-fixed T1 and mask must already be RAS-oriented")
    if t1.shape != mask.shape or not np.allclose(t1.affine, mask.affine, rtol=0, atol=1e-5):
        raise ValueError("held-fixed brain mask must match the T1 grid")
    mask_values = np.asarray(mask.dataobj)
    if not np.isin(mask_values, (0, 1, 255)).all() or not np.any(mask_values):
        raise ValueError("brain mask must be nonempty binary data")
    t1_values = np.asarray(t1.dataobj, dtype=np.float32)
    if not np.isfinite(t1_values).all():
        raise ValueError("T1 contains nonfinite values")
    prepared = np.where(mask_values > 0, t1_values, np.float32(0)).astype("<f4")
    args.raw_output.parent.mkdir(parents=True, exist_ok=True)
    args.raw_output.write_bytes(prepared.tobytes(order="F"))
    result = {
        "scope": "held-fixed input serialization for selected Java/Python TOADS comparison",
        "spatial_convention": "RAS",
        "shape_xyz": list(t1.shape),
        "spacing_mm": [float(v) for v in t1.header.get_zooms()[:3]],
        "array_order": "MIPAV x-fastest (Fortran order), little-endian float32",
        "mask_application": "T1 values outside supplied mask set to zero; same original T1 and mask are passed to Python TOADS",
        "t1_n4": {"path": str(t1_path), "sha256": sha256(t1_path)},
        "brain_mask": {"path": str(mask_path), "sha256": sha256(mask_path)},
        "masked_voxels": int(np.count_nonzero(mask_values)),
        "raw_input": {"path": str(args.raw_output.resolve()), "sha256": sha256(args.raw_output)},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
