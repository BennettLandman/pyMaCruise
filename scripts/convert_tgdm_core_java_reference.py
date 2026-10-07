#!/usr/bin/env python3
"""Restore the synthetic RAS grid metadata on selected Java TGDM outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    args = parser.parse_args()
    grid = nib.load(str(args.input_root / "initial_membership.nii.gz"))
    if nib.aff2axcodes(grid.affine) != ("R", "A", "S"):
        raise ValueError("input grid must be RAS")
    for path in sorted(args.reference_root.glob("*.nii.gz")):
        image = nib.load(str(path))
        if image.shape != grid.shape:
            raise ValueError(f"Java output has wrong shape: {path.name}")
        restored = nib.Nifti1Image(np.asarray(image.dataobj), grid.affine, grid.header.copy())
        nib.save(restored, str(path))
    print(f"Restored RAS geometry on Java outputs in {args.reference_root}")


if __name__ == "__main__":
    main()
