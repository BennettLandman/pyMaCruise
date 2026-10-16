#!/usr/bin/env python3
"""Apply the packet's T1-grid binary brain mask to an N4-corrected NIfTI."""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not args.image.is_file() or not args.mask.is_file():
        parser.error("image and mask must be existing regular files")
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    image = nib.load(str(args.image))
    mask_image = nib.load(str(args.mask))
    if len(image.shape) != 3 or len(mask_image.shape) != 3:
        parser.error("image and mask must both be 3D NIfTI files")
    if image.shape != mask_image.shape or not np.allclose(
        image.affine, mask_image.affine, rtol=0, atol=1e-4
    ):
        parser.error("mask shape and affine must match the N4 image")
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        parser.error("the staged source-run T1 and mask must already be RAS-oriented")

    mask = np.asanyarray(mask_image.dataobj)
    if not np.isfinite(mask).all() or not np.isin(mask, (0, 1)).all():
        parser.error("brain mask must contain finite binary values 0 and 1")
    if not np.any(mask):
        parser.error("brain mask must contain foreground voxels")
    values = np.asanyarray(image.dataobj)
    if not np.isfinite(values).all():
        parser.error("N4 image contains nonfinite voxel values")
    masked = np.where(mask != 0, values, 0).astype(np.float32, copy=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = image.header.copy()
    header.set_data_dtype(np.float32)
    output = nib.Nifti1Image(masked, image.affine, header=header)
    _, qcode = image.get_qform(coded=True)
    _, scode = image.get_sform(coded=True)
    if qcode:
        output.set_qform(image.affine, code=int(qcode), strip_shears=True)
    else:
        output.set_qform(None, code=0)
    output.set_sform(image.affine, code=int(scode) if scode else 1)
    nib.save(output, str(args.output))
    check = nib.load(str(args.output))
    if check.shape != image.shape or not np.allclose(
        check.affine, image.affine, rtol=0, atol=1e-5
    ):
        raise RuntimeError("masked N4 output failed geometry read-back validation")
    print(
        f"Wrote masked N4 image {args.output} shape={check.shape} "
        f"axis_codes={nib.aff2axcodes(check.affine)} "
        f"masked_voxels={int(np.count_nonzero(mask == 0))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
