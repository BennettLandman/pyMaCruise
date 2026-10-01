from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent


def main() -> None:
    size = 41
    axis = np.arange(size, dtype=np.float32) - np.float32(size // 2)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.sqrt(x * x + y * y + z * z)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, np.hypot(x, y))

    # A closed, cortex-like star-shaped phantom with intersecting angular
    # frequencies that create repeated sulcal/gyral folds and saddle cells.
    fold = (
        np.float32(1.15) * np.cos(np.float32(5.0) * theta) * np.cos(np.float32(3.0) * phi)
        + np.float32(0.45) * np.sin(np.float32(8.0) * theta) * np.sin(np.float32(4.0) * phi)
    )
    base = (radius - (np.float32(12.0) + fold)).astype(np.float32)
    fields = {
        "inner": base + np.float32(2.0),
        "central": base,
        "outer": base - np.float32(2.0),
    }
    affine = np.eye(4, dtype=np.float64)
    for name, values in fields.items():
        image = nib.Nifti1Image(values, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, HERE / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
