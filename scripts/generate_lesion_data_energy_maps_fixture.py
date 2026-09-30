#!/usr/bin/env python3
"""Generate the fixed multi-channel RAS volume used by P39zzzzm tests."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/lesion_data_energy_maps"


def main() -> None:
    FIXTURE.mkdir(parents=True, exist_ok=True)
    nx, ny, nz, channels = 5, 6, 7, 3
    x, y, z, c = np.indices((nx, ny, nz, channels), dtype=np.int32)
    image = ((13 * x - 7 * y + 5 * z + 11 * c - 23).astype(np.float32) / np.float32(8.0))
    affine = np.diag([1.5, 2.0, 2.5, 1.0]).astype(np.float64)
    nii = nib.Nifti1Image(image, affine)
    nii.set_qform(affine, code=1)
    nii.set_sform(affine, code=1)
    nib.save(nii, FIXTURE / "images_ras.nii.gz")


if __name__ == "__main__":
    main()
