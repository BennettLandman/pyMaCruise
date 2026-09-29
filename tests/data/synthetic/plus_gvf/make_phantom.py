from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent


def main() -> None:
    axis = np.linspace(-1.0, 1.0, 13, dtype=np.float32)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    values = (
        30.0
        + 110.0 * np.exp(-((x / 0.62) ** 2 + (y / 0.35) ** 2 + (z / 0.48) ** 2))
        + 22.0 * np.exp(-(((x - 0.38) / 0.18) ** 2 + ((y + 0.2) / 0.25) ** 2 + (z / 0.3) ** 2))
    ).astype(np.float32)
    image = nib.Nifti1Image(values, np.diag([1.0, 1.2, 1.5, 1.0]))
    image.set_qform(image.affine, code=1)
    image.set_sform(image.affine, code=1)
    nib.save(image, HERE / "scalar.nii.gz")


if __name__ == "__main__":
    main()
