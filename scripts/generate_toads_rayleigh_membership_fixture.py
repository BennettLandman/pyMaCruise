"""Generate the folded RAS inputs/maps for the RAYLEIGH pass composition."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_rayleigh_membership"


def _write(name: str, values: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / name)


def main() -> None:
    nx, ny, nz, channels, classes = 9, 8, 7, 3, 3
    x, y, z = np.ogrid[:nx, :ny, :nz]
    dx = (x - 4.2) / 1.25
    dy = (y - 3.15) / 1.05
    dz = (z - 2.7) / 1.0
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    angle = np.arctan2(dy, dx)
    inner = 1.35 + 0.27 * np.sin(4 * angle + dz * 0.7)
    outer = 3.25 + 0.48 * np.sin(5 * angle + dz * 0.4 + 0.35)
    ribbon = ((radius >= inner) & (radius <= outer)).astype(np.uint8)

    images = np.empty((nx, ny, nz, channels), dtype=np.float32)
    field = np.empty_like(images)
    bias = np.empty_like(images)
    for c in range(channels):
        code = (5 * x + 7 * y + 11 * z + 3 * c) % 13
        images[..., c] = np.asarray(
            np.float32(-0.18 + 0.045 * c)
            + np.asarray(x, dtype=np.float32) * np.float32(0.038)
            + np.asarray(y, dtype=np.float32) * np.float32(0.046)
            + np.asarray(z, dtype=np.float32) * np.float32(0.052)
            + np.asarray(ribbon, dtype=np.float32) * np.float32(0.09)
            + np.asarray(code, dtype=np.float32) * np.float32(0.006),
            dtype=np.float32,
        )
        field_code = (3 * x + 4 * y + 5 * z + c) % 7
        field[..., c] = np.asarray(
            np.float32(0.82) + np.asarray(field_code, dtype=np.float32) * np.float32(0.035),
            dtype=np.float32,
        )
        bias_code = (x + 2 * y + 3 * z + 2 * c) % 9
        bias[..., c] = np.asarray(
            np.float32(-0.025) + np.asarray(bias_code, dtype=np.float32) * np.float32(0.012),
            dtype=np.float32,
        )

    memberships = np.empty((nx, ny, nz, classes), dtype=np.float32)
    shapes = np.empty_like(memberships)
    labels = np.full((nx, ny, nz, 4), -1, dtype=np.int16)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for k in range(classes):
                    mi = (ix + 2 * iy + iz + 3 * k) % 7
                    memberships[ix, iy, iz, k] = np.float32(
                        np.float32(0.14) + np.float32(np.float32(mi) * np.float32(0.019))
                    )
                    si = (2 * ix + iy + 3 * iz + k) % 5
                    shapes[ix, iy, iz, k] = np.float32(
                        np.float32(0.22)
                        + np.float32(np.float32(0.47) * ribbon[ix, iy, iz])
                        + np.float32(np.float32(si) * np.float32(0.013))
                    )
                labels[ix, iy, iz, 0] = (ix + iy + iz) % classes
                if (ix + 2 * iy + iz) % 3 == 0:
                    labels[ix, iy, iz, 1] = (ix + iy + iz + 1) % classes
                if (ix + iy + 2 * iz) % 2 == 0:
                    labels[ix, iy, iz, 2] = (ix + 2 * iy + iz + 1) % classes
                if (ix + 2 * iy + 2 * iz) % 3 == 0:
                    labels[ix, iy, iz, 3] = (ix + 2 * iy + iz + 2) % classes

    affine = np.asarray(
        [[0.9, 0, 0, -4.05], [0, 1.2, 0, -4.2], [0, 0, 1.45, -4.35], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, array in (
        ("images.nii.gz", images),
        ("field.nii.gz", field),
        ("psf_bias.nii.gz", bias),
        ("folded_ribbon.nii.gz", ribbon),
        ("initial_memberships.nii.gz", memberships),
        ("transformed_shapes.nii.gz", shapes),
        ("closest_labels.nii.gz", labels),
    ):
        _write(name, array, affine)
    print(f"wrote folded RAS membership-pass inputs to {OUTPUT}")


if __name__ == "__main__":
    main()
