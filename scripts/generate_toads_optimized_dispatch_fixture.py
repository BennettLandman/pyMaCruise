"""Generate the deterministic mixed optimized-tissue TOADS RAS fixture."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_optimized_dispatch"


def _write(name: str, values: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / name)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nx, ny, nz, channels = 15, 13, 11, 3
    x, y, z = np.ogrid[:nx, :ny, :nz]
    dx = (x - 7.2) / 1.05
    dy = (y - 5.35) / 0.95
    dz = (z - 4.45) / 1.1
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    angle = np.arctan2(dy, dx)
    inner = 3.2 + 0.65 * np.sin(5 * angle + 0.45 * dz)
    outer = 5.55 + 0.72 * np.sin(5 * angle + 0.35 * dz + 0.4)
    ribbon = ((radius >= inner) & (radius <= outer)).astype(np.uint8)

    image = np.empty((nx, ny, nz, channels), dtype=np.float32)
    field = np.empty_like(image)
    bias = np.empty_like(image)
    for c in range(channels):
        code = (7 * x + 11 * y + 13 * z + 5 * c) % 29
        tissue = np.float32(0.18 + 0.014 * c) + np.float32(0.46) * ribbon
        variation = np.asarray(code, dtype=np.float32) * np.float32(0.007)
        image[..., c] = np.asarray(tissue + variation, dtype=np.float32)
        field_code = (3 * x + 5 * y + 7 * z + c) % 9
        field[..., c] = np.asarray(np.float32(0.82) + np.asarray(field_code, dtype=np.float32) * np.float32(0.025), dtype=np.float32)
        bias_code = (x + 2 * y + 3 * z + 2 * c) % 11
        bias[..., c] = np.asarray(np.asarray(bias_code, dtype=np.float32) * np.float32(0.003), dtype=np.float32)

    affine = np.asarray(
        [[1.05, 0.0, 0.0, -7.35], [0.0, 1.2, 0.0, -7.2], [0.0, 0.0, 1.35, -6.75], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    _write("images.nii.gz", image, affine)
    _write("field.nii.gz", field, affine)
    _write("psf_bias.nii.gz", bias, affine)
    _write("folded_ribbon.nii.gz", ribbon, affine)
    print(f"wrote fixture to {OUTPUT}")


if __name__ == "__main__":
    main()
