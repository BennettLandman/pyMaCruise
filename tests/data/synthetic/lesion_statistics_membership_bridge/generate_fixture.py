"""Generate deterministic folded inputs for the P39zzzzs bridge harness."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (6, 5, 4)
CHANNELS = 2
CLASSES = 4


def make_inputs() -> np.ndarray:
    """Return RAS frames: images, fields, labels, memberships, auxiliary, shapes."""
    nx, ny, nz = SHAPE
    values = np.zeros((*SHAPE, 16), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                fold = (x * 11 + y * 7 + z * 3 + x * y + 2 * y * z) % 23
                for c in range(CHANNELS):
                    values[x, y, z, c] = np.float32(
                        np.float32(0.15) + np.float32((fold + 5 * c) % 23) * np.float32(0.025)
                    )
                    values[x, y, z, 2 + c] = np.float32(
                        np.float32(0.82) + np.float32((x * 2 + y * 3 + z * 4 + c) % 5) * np.float32(0.04)
                    )
                values[x, y, z, 4] = np.float32((x + 2 * y + z + (x * y + 2 * z) % 3) % CLASSES)
                for k in range(CLASSES):
                    values[x, y, z, 5 + k] = np.float32(
                        np.float32(0.08) + np.float32((fold + 3 * k) % 9) * np.float32(0.045)
                    )
                    values[x, y, z, 12 + k] = np.float32(
                        np.float32(0.45) + np.float32((x * 3 + y * 5 + z * 7 + k * 2) % 11) * np.float32(0.035)
                    )
                values[x, y, z, 9] = np.float32(
                    np.float32(0.04) + np.float32((x + 2 * y + z) % 8) * np.float32(0.06)
                )
                values[x, y, z, 10] = np.float32(
                    np.float32(0.03) + np.float32((2 * x + y + 3 * z) % 8) * np.float32(0.065)
                )
                values[x, y, z, 11] = np.float32((x + 3 * y + z) % 5) * np.float32(0.12)
    affine = np.eye(4, dtype=np.float64)
    affine[:3, 3] = (-17.0, 9.0, 4.5)
    image = nib.Nifti1Image(values, affine)
    image.set_data_dtype(np.float32)
    nib.save(image, ROOT / "bridge_inputs_ras.nii.gz")


if __name__ == "__main__":
    make_inputs()
