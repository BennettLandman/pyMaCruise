from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (15, 15, 15)
CLASSES = 4


def main() -> None:
    values = np.zeros((*SHAPE, 19), dtype=np.float32)
    values[..., 17] = np.float32(-1.0)
    nx, ny, nz = SHAPE
    segmentation = np.empty(SHAPE, dtype=np.int16)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                fold = (x * 11 + y * 7 + z * 3 + x * y + 2 * y * z) % 23
                warp = ((3 * y + 2 * z + y * z) % 5) - 2
                if x < 4 + warp:
                    label = 0
                elif x < 8 + warp:
                    label = 1
                elif x < 12 + warp:
                    label = 2
                else:
                    label = 3
                segmentation[x, y, z] = label
                for channel in range(2):
                    values[x, y, z, channel] = np.float32(
                        np.float32(0.12 + 0.19 * label)
                        + np.float32((fold + 5 * channel) % 17) * np.float32(0.011)
                    )
                    values[x, y, z, 2 + channel] = np.float32(
                        0.81 + 0.035 * ((2 * x + 3 * y + 5 * z + channel) % 5)
                    )
                values[x, y, z, 4] = np.float32(label)
                for k in range(CLASSES):
                    values[x, y, z, 5 + k] = np.float32(
                        0.12 + 0.037 * ((fold + 3 * k + x * y) % 13)
                    )
                    values[x, y, z, 12 + k] = np.float32(
                        0.42 + 0.027 * ((x * 3 + y * 5 + z * 7 + 2 * k) % 15)
                    )
                values[x, y, z, 9] = np.float32(0.08 + 0.055 * ((x + 2 * y + z) % 8))
                values[x, y, z, 10] = np.float32(0.06 + 0.05 * ((2 * x + y + 3 * z) % 9))
                values[x, y, z, 11] = np.float32(0.03 + 0.10 * ((x + 3 * y + z) % 6))
                values[x, y, z, 16] = np.float32(label)
                values[x, y, z, 18] = np.float32(-((13 * x + 7 * y + 3 * z) % 19) / 10.0)
    values[4, 4, 4, 17] = values[4, 4, 4, 16]
    values[2, 5, 6, 17] = values[2, 5, 6, 16]
    affine = np.eye(4, dtype=np.float64)
    affine[:3, 3] = (-18.0, 6.5, 11.0)
    image = nib.Nifti1Image(values, affine)
    image.set_data_dtype(np.float32)
    nib.save(image, str(ROOT / "iteration_inputs_ras.nii.gz"))


if __name__ == "__main__":
    main()
