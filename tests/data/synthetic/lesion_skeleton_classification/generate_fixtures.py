from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def _save(name: str, values: list[int]) -> None:
    array = np.asarray(values, dtype=np.int8).reshape((10, 1, 1))
    image = nib.Nifti1Image(array, AFFINE)
    image.set_qform(AFFINE, code=1)
    image.set_sform(AFFINE, code=1)
    nib.save(image, HERE / f"{name}.nii.gz")


def main() -> None:
    _save("initial_classification", [5, -10, 20, 0, 99, 0, 2, -10, -7, 20])
    _save("skeleton", [-1, 0, 1, 2, 3, -1, 2, 0, 3, 1])


if __name__ == "__main__":
    main()
