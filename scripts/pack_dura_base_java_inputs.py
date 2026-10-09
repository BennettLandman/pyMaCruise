"""Pack four RAS NIfTI images for the selected Java base-Dura harness."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import nibabel as nib
import numpy as np


def pack(inputs: list[Path], output: Path) -> tuple[int, int, int]:
    reference = nib.load(str(inputs[0]))
    shape = tuple(int(value) for value in reference.shape)
    if len(shape) != 3:
        raise ValueError("Dura Java harness requires three-dimensional images")
    if nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError(f"T1 input is not RAS-oriented: {inputs[0]}")
    arrays: list[np.ndarray] = []
    for path in inputs:
        image = nib.load(str(path))
        if image.shape != shape or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"input is not on the common RAS grid: {path}")
        if not np.allclose(image.affine, reference.affine, atol=1e-5):
            raise ValueError(f"input affine does not match T1: {path}")
        values = np.asarray(image.dataobj, dtype=np.float32)
        if not np.isfinite(values).all():
            raise ValueError(f"input contains non-finite values: {path}")
        arrays.append(values.astype(">f4", copy=False))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(
        struct.pack(">III", *shape)
        + b"".join(values.tobytes(order="C") for values in arrays)
    )
    return shape


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("t1", type=Path)
    parser.add_argument("white_matter", type=Path)
    parser.add_argument("gray_matter", type=Path)
    parser.add_argument("white_matter_mask", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    shape = pack(
        [args.t1, args.white_matter, args.gray_matter, args.white_matter_mask],
        args.output,
    )
    print(f"wrote {args.output} for {shape} RAS voxels")


if __name__ == "__main__":
    main()
