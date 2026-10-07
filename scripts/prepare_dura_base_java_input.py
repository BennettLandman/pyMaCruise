"""Pack held-fixed brain-resolution RAS NIfTIs for the Java Dura harness."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import nibabel as nib
import numpy as np


NAMES = ("t1", "white_matter", "gray_matter", "white_matter_mask")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    images = [nib.load(args.input_dir / f"{name}.nii.gz") for name in NAMES]
    shape = images[0].shape
    affine = images[0].affine
    if len(shape) != 3 or any(image.shape != shape for image in images):
        raise ValueError("Dura Java inputs must share one 3D grid")
    if nib.aff2axcodes(affine) != ("R", "A", "S"):
        raise ValueError("Dura Java inputs must be RAS-oriented")
    if any(not np.array_equal(image.affine, affine) for image in images[1:]):
        raise ValueError("Dura Java inputs must have identical affines")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for image in images:
            values = np.asarray(image.dataobj, dtype=np.float32)
            if not np.isfinite(values).all():
                raise ValueError("Dura Java inputs must be finite")
            values.astype(">f4", copy=False).tofile(stream)
    print(f"packed {len(images)} RAS float32 volumes with shape {shape} to {args.output}")


if __name__ == "__main__":
    main()
