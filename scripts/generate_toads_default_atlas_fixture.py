"""Generate the selected default-atlas 11-class synthetic RAS fixture."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_default_atlas_dispatch"
CLASSES = 11


def _write(name: str, data: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(data, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / name)


def main() -> None:
    nx, ny, nz, channels = 9, 8, 7, 2
    x, y, z = np.ogrid[:nx, :ny, :nz]
    dx, dy, dz = (x - 4.1) / 1.3, (y - 3.2) / 1.05, (z - 2.65) / .95
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    angle = np.arctan2(dy, dx)
    inner = 1.25 + .24 * np.sin(4 * angle + .65 * dz)
    outer = 3.15 + .44 * np.sin(5 * angle + .4 * dz + .3)
    ribbon = ((radius >= inner) & (radius <= outer)).astype(np.uint8)
    images = np.empty((nx, ny, nz, channels), dtype=np.float32)
    field = np.empty_like(images)
    bias = np.empty_like(images)
    for c in range(channels):
        code = (5 * x + 7 * y + 11 * z + 3 * c) % 13
        images[..., c] = np.asarray(
            np.float32(-.12 + .035 * c)
            + np.asarray(x, dtype=np.float32) * np.float32(.031)
            + np.asarray(y, dtype=np.float32) * np.float32(.037)
            + np.asarray(z, dtype=np.float32) * np.float32(.044)
            + ribbon.astype(np.float32) * np.float32(.075)
            + np.asarray(code, dtype=np.float32) * np.float32(.005), dtype=np.float32
        )
        fcode = (3 * x + 5 * y + 2 * z + c) % 7
        field[..., c] = np.asarray(
            np.float32(.83) + np.asarray(fcode, dtype=np.float32) * np.float32(.037), dtype=np.float32
        )
        bcode = (x + 2 * y + 3 * z + 2 * c) % 9
        bias[..., c] = np.asarray(
            np.float32(-.025) + np.asarray(bcode, dtype=np.float32) * np.float32(.014), dtype=np.float32
        )

    initial = np.empty((nx, ny, nz, CLASSES), dtype=np.float32)
    shapes = np.empty_like(initial)
    labels = np.full((nx, ny, nz, 4), -1, dtype=np.int16)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for k in range(CLASSES):
                    mi = (ix + 2 * iy + iz + 3 * k) % 7
                    initial[ix, iy, iz, k] = np.float32(
                        np.float32(.13) + np.float32(np.float32(mi) * np.float32(.017))
                    )
                    si = (2 * ix + iy + 3 * iz + k) % 5
                    shapes[ix, iy, iz, k] = np.float32(
                        np.float32(.24) + np.float32(np.float32(.43) * ribbon[ix, iy, iz])
                        + np.float32(np.float32(si) * np.float32(.012))
                    )
                labels[ix, iy, iz, 0] = (ix + iy + iz) % CLASSES
                if (ix + 2 * iy + iz) % 3 == 0:
                    labels[ix, iy, iz, 1] = (labels[ix, iy, iz, 0] + 2) % CLASSES
                if (ix + iy + 2 * iz) % 2 == 0:
                    labels[ix, iy, iz, 2] = (labels[ix, iy, iz, 0] + 1) % CLASSES
                if (ix + 2 * iy + 2 * iz) % 3 == 0:
                    labels[ix, iy, iz, 3] = (labels[ix, iy, iz, 0] + 3) % CLASSES

    affine = np.asarray(
        [[.9, 0, 0, -4.05], [0, 1.2, 0, -4.2], [0, 0, 1.45, -4.35], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, data in (
        ("images.nii.gz", images), ("field.nii.gz", field), ("psf_bias.nii.gz", bias),
        ("folded_ribbon.nii.gz", ribbon), ("initial_memberships.nii.gz", initial),
        ("transformed_shapes.nii.gz", shapes), ("closest_labels.nii.gz", labels),
    ):
        _write(name, data, affine)
    print(f"wrote folded true-RAS 11-class fixture to {OUTPUT}")


if __name__ == "__main__":
    main()
