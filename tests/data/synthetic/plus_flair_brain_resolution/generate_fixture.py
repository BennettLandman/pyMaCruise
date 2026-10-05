from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import shift


ROOT = Path(__file__).parent
SHAPE = (172, 216, 164)


def main() -> None:
    x = (np.arange(SHAPE[0], dtype=np.float32) - np.float32((SHAPE[0] - 1) / 2))[:, None, None]
    y = (np.arange(SHAPE[1], dtype=np.float32) - np.float32((SHAPE[1] - 1) / 2))[None, :, None]
    z = (np.arange(SHAPE[2], dtype=np.float32) - np.float32((SHAPE[2] - 1) / 2))[None, None, :]
    radius = np.sqrt((x / 67.0) ** 2 + (y / 86.0) ** 2 + (z / 65.0) ** 2)
    angle = np.arctan2(y / 86.0, x / 67.0)
    folded_boundary = 0.78 + 0.025 * np.sin(9.0 * angle + z / 12.0)
    brain = radius < folded_boundary
    deep = (x / 29.0) ** 2 + (y / 39.0) ** 2 + (z / 32.0) ** 2 < 1.0
    ventricles = (x / 10.0) ** 2 + (y / 21.0) ** 2 + (z / 13.0) ** 2 < 1.0

    t1 = np.zeros(SHAPE, dtype=np.float32)
    t1[brain] = np.float32(0.23)
    t1[deep & brain] = np.float32(0.62)
    t1[ventricles] = np.float32(0.12)
    t1 += (0.025 * np.cos(0.11 * x) * np.cos(0.085 * y) * np.cos(0.13 * z)).astype(np.float32)
    t1 = np.clip(t1, 0.0, 1.0).astype(np.float32, copy=False)

    flair = np.sqrt(t1, dtype=np.float32)
    focal = ((x - 18.0) / 8.0) ** 2 + ((y + 24.0) / 11.0) ** 2 + ((z - 9.0) / 7.0) ** 2 < 1.0
    flair[focal] = np.float32(1.0)
    # The moving image is displaced by (2, -3, 1) voxels; registration maps it
    # back to the fixed T1 grid.
    flair = shift(flair, shift=(2.0, -3.0, 1.0), order=1, mode="constant", cval=0.0, prefilter=False)

    affine = np.eye(4, dtype=np.float64)
    for name, values in (("t1.nii.gz", t1), ("flair_shifted.nii.gz", flair)):
        image = nib.Nifti1Image(values.astype(np.float32, copy=False), affine)
        image.set_data_dtype(np.float32)
        nib.save(image, str(ROOT / name))


if __name__ == "__main__":
    main()
