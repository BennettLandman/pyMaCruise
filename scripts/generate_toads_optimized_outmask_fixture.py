"""Generate folded asymmetric true-RAS inputs for P18cze."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_optimized_outmask"
TYPES = ("OPTIMIZED_GM", "OUTMASK", "MASK", "OBJECT", "OPTIMIZED_CSF", "OPTIMIZED_WM")


def _save(name: str, data: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(data, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / name)


def main() -> None:
    nx, ny, nz, channels, classes = 9, 8, 7, 2, len(TYPES)
    x, y, z = np.ogrid[:nx, :ny, :nz]
    dx, dy, dz = (x - 4.1) / 1.3, (y - 3.2) / 1.05, (z - 2.65) / 0.95
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    angle = np.arctan2(dy, dx)
    inner = 1.25 + 0.24 * np.sin(4 * angle + 0.65 * dz)
    outer = 3.15 + 0.44 * np.sin(5 * angle + 0.4 * dz + 0.3)
    ribbon = ((radius >= inner) & (radius <= outer)).astype(np.uint8)
    images = np.empty((nx, ny, nz, channels), dtype=np.float32)
    field = np.empty_like(images)
    bias = np.empty_like(images)
    for c in range(channels):
        code = (5 * x + 7 * y + 11 * z + 3 * c) % 13
        images[..., c] = np.asarray(
            np.float32(-0.12 + 0.035 * c)
            + np.asarray(x, dtype=np.float32) * np.float32(0.031)
            + np.asarray(y, dtype=np.float32) * np.float32(0.037)
            + np.asarray(z, dtype=np.float32) * np.float32(0.044)
            + ribbon.astype(np.float32) * np.float32(0.075)
            + np.asarray(code, dtype=np.float32) * np.float32(0.005),
            dtype=np.float32,
        )
        fcode = (3 * x + 5 * y + 2 * z + c) % 7
        field[..., c] = np.asarray(
            np.float32(0.83) + np.asarray(fcode, dtype=np.float32) * np.float32(0.037),
            dtype=np.float32,
        )
        bcode = (x + 2 * y + 3 * z + 2 * c) % 9
        bias[..., c] = np.asarray(
            np.float32(-0.025) + np.asarray(bcode, dtype=np.float32) * np.float32(0.014),
            dtype=np.float32,
        )

    initial = np.empty((nx, ny, nz, classes), dtype=np.float32)
    shapes = np.empty_like(initial)
    labels = np.full((nx, ny, nz, 4), -1, dtype=np.int16)
    outlier = np.empty((nx, ny, nz), dtype=np.float32)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for k in range(classes):
                    mi = (ix + 2 * iy + iz + 3 * k) % 7
                    initial[ix, iy, iz, k] = np.float32(
                        np.float32(0.13) + np.float32(np.float32(mi) * np.float32(0.017))
                    )
                    si = (2 * ix + iy + 3 * iz + k) % 5
                    shapes[ix, iy, iz, k] = np.float32(
                        np.float32(0.24)
                        + np.float32(np.float32(0.43) * ribbon[ix, iy, iz])
                        + np.float32(np.float32(si) * np.float32(0.012))
                    )
                labels[ix, iy, iz, 0] = (ix + iy + iz) % classes
                if (ix + 2 * iy + iz) % 3 == 0:
                    labels[ix, iy, iz, 1] = (ix + iy + iz + 1) % classes
                if (ix + iy + 2 * iz) % 2 == 0:
                    labels[ix, iy, iz, 2] = (ix + 2 * iy + iz + 1) % classes
                if (ix + 2 * iy + 2 * iz) % 3 == 0:
                    labels[ix, iy, iz, 3] = (ix + 2 * iy + iz + 2) % classes
                oi = (2 * ix + 3 * iy + iz) % 7
                outlier[ix, iy, iz] = np.float32(
                    np.float32(0.18) + np.float32(np.float32(oi) * np.float32(0.055))
                )

    affine = np.asarray(
        [[0.9, 0, 0, -4.05], [0, 1.2, 0, -4.2], [0, 0, 1.45, -4.35], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, data in (
        ("images.nii.gz", images), ("field.nii.gz", field), ("psf_bias.nii.gz", bias),
        ("folded_ribbon.nii.gz", ribbon), ("initial_memberships.nii.gz", initial),
        ("transformed_shapes.nii.gz", shapes), ("closest_labels.nii.gz", labels),
        ("initial_outlier_membership.nii.gz", outlier),
    ):
        _save(name, data, affine)
    print(f"wrote folded true-RAS P18cze fixture to {OUTPUT}")


if __name__ == "__main__":
    main()
