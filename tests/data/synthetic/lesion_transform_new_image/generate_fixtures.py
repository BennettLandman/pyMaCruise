from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent


def main() -> None:
    x, y, z = np.indices((9, 10, 11), dtype=np.int32)
    values = ((x * 7 + y * 3 + z * 5 + x * y + z * y * 2) % 97).astype(np.float32)
    values = np.asarray(values / np.float32(16.0), dtype=np.float32)
    affine = np.asarray(
        [[0.9, 0, 0, 12], [0, 1.0, 0, -7], [0, 0, 1.1, 4], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, HERE / "template_image.nii.gz")


if __name__ == "__main__":
    main()
