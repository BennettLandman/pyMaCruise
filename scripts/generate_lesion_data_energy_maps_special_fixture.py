#!/usr/bin/env python3
"""Generate folded multi-channel RAS inputs for P39zzzzn."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/lesion_data_energy_maps_special"


def main() -> None:
    FIXTURE.mkdir(parents=True, exist_ok=True)
    nx, ny, nz, channels = 7, 8, 9, 3
    x, y, z, c = np.indices((nx, ny, nz, channels), dtype=np.int32)
    fold = (3 * x + 2 * y + z) % 5 - 2
    intensity = (13 * x - 9 * y + 7 * z + 5 * fold + 11 * c - 20).astype(np.float32) / np.float32(32.0)
    field = np.float32(0.625) + ((x + 2 * y + 3 * z + c) % 7).astype(np.float32) / np.float32(16.0)
    combined = np.concatenate((intensity, field), axis=-1)
    affine = np.diag([1.5, 2.0, 2.5, 1.0]).astype(np.float64)
    nii = nib.Nifti1Image(combined, affine)
    nii.set_qform(affine, code=1)
    nii.set_sform(affine, code=1)
    nib.save(nii, FIXTURE / "inputs_ras.nii.gz")


if __name__ == "__main__":
    main()
