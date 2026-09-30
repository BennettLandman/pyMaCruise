"""Create deterministic active-class inputs for the P39zzzzp Java harness."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (5, 4, 5)
CLASSES = 4


def _f32(value: int, denominator: int, offset_num: int = 0, offset_den: int = 1) -> np.float32:
    return np.float32(np.float32(value) / np.float32(denominator) + np.float32(offset_num) / np.float32(offset_den))


def make_inputs() -> np.ndarray:
    """Return ``(x,y,z,18)`` RAS float32 inputs in documented channel order."""
    xsize, ysize, zsize = SHAPE
    values = np.zeros((*SHAPE, 18), dtype=np.float32)
    for x in range(xsize):
        for y in range(ysize):
            for z in range(zsize):
                # The asymmetric modular folds make neighboring energy and
                # shape values vary independently on all three axes.
                fold = (x * 7 + y * 11 + z * 5 + x * y) % 29
                for k in range(CLASSES):
                    values[x, y, z, k] = _f32((fold + 3 * k) % 17 + 1, 32)
                    values[x, y, z, 7 + k] = _f32((fold * 3 + 5 * k) % 23 + 4, 16)
                    values[x, y, z, 14 + k] = _f32((x * 3 + y * 5 + z * 7 + k * 2) % 12 + 4, 16)
                values[x, y, z, 4] = _f32((fold * 5) % 13 + 2, 16)
                values[x, y, z, 5] = _f32((fold * 7) % 11 + 2, 16)
                values[x, y, z, 6] = _f32((fold * 2) % 7 + 1, 16)
                values[x, y, z, 11] = _f32((fold * 9) % 19 + 3, 16)
                values[x, y, z, 12] = _f32((fold * 4) % 17 + 2, 16)
                values[x, y, z, 13] = _f32((fold * 6) % 15 + 2, 16)
    return values


def main() -> None:
    affine = np.eye(4, dtype=np.float64)
    affine[:3, 3] = (-23.0, 14.0, 8.0)
    image = nib.Nifti1Image(make_inputs(), affine)
    image.set_data_dtype(np.float32)
    nib.save(image, ROOT / "inputs_ras.nii.gz")


if __name__ == "__main__":
    main()
