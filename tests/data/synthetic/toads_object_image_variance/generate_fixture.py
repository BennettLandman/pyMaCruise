from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (4, 3, 3)
AFFINE = np.eye(4, dtype=np.float64)


def main() -> None:
    images = np.empty(SHAPE + (2,), dtype=np.float32)
    field = np.ones_like(images)
    memberships = np.empty(SHAPE + (3,), dtype=np.float32)
    for x in range(SHAPE[0]):
        for y in range(SHAPE[1]):
            for z in range(SHAPE[2]):
                images[x, y, z, 0] = np.float32(0.5)
                images[x, y, z, 1] = np.float32(
                    0.12 + ((7 * x + 5 * y + 3 * z) % 19) * 0.027
                )
                field[x, y, z, 1] = np.float32(
                    0.8 + ((x + 3 * y + 5 * z) % 7) * 0.04
                )
                for k in range(3):
                    memberships[x, y, z, k] = np.float32(
                        ((3 * x + y + 2 * z + 2 * k) % 7) / 6.0
                    )
    for name, array in (
        ("images", images),
        ("memberships", memberships),
        ("intensity_field", field),
    ):
        image = nib.Nifti1Image(array, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
