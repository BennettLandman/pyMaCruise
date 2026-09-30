from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent


def main() -> None:
    ordinary = np.asarray(
        [
            [0.2, 0.4, 0.3], [0.5, 0.5, 0.1], [0, 0, 0], [-1, -2, -3],
            [np.nan, 0.1, 0], [np.inf, np.inf, 0], [0.8, 0.2, 0],
            [0, 0, 0], [0, 0.4, 0.4], [0, 0, 0],
        ],
        dtype=np.float32,
    ).reshape((10, 1, 1, 3))
    lesion = np.asarray([0.5, 0.5, 0, -0.5, 0.2, np.inf, 0.8, 0.01, 0.1, np.nan], dtype=np.float32)
    lesion = lesion.reshape((10, 1, 1, 1))
    values = np.concatenate((ordinary, lesion), axis=-1)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, HERE / "memberships_and_lesion.nii.gz")


if __name__ == "__main__":
    main()
