"""Generate asymmetric true-RAS NIfTI inputs for the RAYLEIGH energy port."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_rayleigh_energy"


def _write(name: str, values: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / name)


def main() -> None:
    nx, ny, nz, channels = 11, 9, 7, 3
    x, y, z = np.ogrid[:nx, :ny, :nz]
    image = np.empty((nx, ny, nz, channels), dtype=np.float32)
    field = np.empty_like(image)
    bias = np.empty_like(image)
    for c in range(channels):
        pattern = (5 * x + 7 * y + 11 * z + 3 * c) % 13
        image[..., c] = np.asarray(
            np.float32(-0.28 + 0.065 * c)
            + np.asarray(x, dtype=np.float32) * np.float32(0.041)
            + np.asarray(y, dtype=np.float32) * np.float32(0.053)
            + np.asarray(z, dtype=np.float32) * np.float32(0.047)
            + np.asarray(pattern, dtype=np.float32) * np.float32(0.009),
            dtype=np.float32,
        )
        field_code = (3 * x + 5 * y + 2 * z + c) % 7
        field[..., c] = np.asarray(
            np.float32(0.78) + np.asarray(field_code, dtype=np.float32) * np.float32(0.055),
            dtype=np.float32,
        )
        bias_code = (x + 2 * y + 4 * z + 3 * c) % 9
        bias[..., c] = np.asarray(
            np.float32(-0.035) + np.asarray(bias_code, dtype=np.float32) * np.float32(0.018),
            dtype=np.float32,
        )

    affine = np.asarray(
        [[0.85, 0, 0, -4.25], [0, 1.15, 0, -4.6], [0, 0, 1.4, -4.2], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    _write("images.nii.gz", image, affine)
    _write("field.nii.gz", field, affine)
    _write("psf_bias.nii.gz", bias, affine)
    print(f"wrote asymmetric RAS energy inputs to {OUTPUT}")


if __name__ == "__main__":
    main()
