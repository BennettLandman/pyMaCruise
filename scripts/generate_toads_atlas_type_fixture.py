"""Generate the committed folded RAS fixture for P18czk."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_atlas_type_membership"
SUBJECT_SHAPE = (9, 8, 7)
TEMPLATE_SHAPE = (11, 10, 9)
CLASSES = 13
CHANNELS = 2
AFFINE = np.asarray(
    [[0.9, 0, 0, -31.25], [0, 1.2, 0, -42.5], [0, 0, 1.45, 18.75], [0, 0, 0, 1]],
    dtype=np.float32,
)


def _write(name: str, values: np.ndarray, affine: np.ndarray = AFFINE) -> None:
    image = nib.Nifti1Image(np.asarray(values), affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, str(ROOT / f"{name}.nii.gz"))


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    nx, ny, nz = SUBJECT_SHAPE
    grid = np.indices(SUBJECT_SHAPE, dtype=np.float32)
    x, y, z = grid
    folded = ((x - 4.0) ** 2 / 10 + (y - 3.0) ** 2 / 6 + (z - 3.0) ** 2 / 5 < 1)
    folded ^= ((x + y + z) % 5 == 0)
    image = np.empty((*SUBJECT_SHAPE, CHANNELS), dtype=np.float32)
    image[..., 0] = 0.17 + 0.031 * x + 0.019 * y + 0.013 * z + 0.07 * folded
    image[..., 1] = 0.11 + 0.014 * x + 0.027 * y + 0.021 * z + 0.09 * (folded & ((x + z) % 2 == 0))
    image[0, 2, 3, 0] = 0.0
    image[-1, 4, 2, 1] = -0.0

    field = np.empty_like(image)
    bias = np.empty_like(image)
    for channel in range(CHANNELS):
        field[..., channel] = 0.91 + 0.006 * x + 0.004 * y + 0.003 * z + 0.01 * channel
        bias[..., channel] = 0.002 * (x - 3) - 0.001 * (y + 1) + 0.0015 * z + .0007 * channel

    initial = np.empty((*SUBJECT_SHAPE, CLASSES), dtype=np.float32)
    shapes = np.empty_like(initial)
    for k in range(CLASSES):
        values = .035 + .002 * ((x + 2 * y + z + k) % 7)
        initial[..., k] = values
        shape = .19 + .012 * ((2 * x + y + 3 * z + k) % 9)
        shapes[..., k] = shape

    closest = np.full((*SUBJECT_SHAPE, 4), -1, dtype=np.int16)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                closest[ix, iy, iz] = (
                    (ix + 2 * iy + iz) % CLASSES,
                    (ix + iy + 3 * iz + 4) % CLASSES,
                    (2 * ix + iy + iz + 1) % CLASSES,
                    (ix + 3 * iy + iz + 7) % CLASSES,
                )

    template_grid = np.indices(TEMPLATE_SHAPE, dtype=np.float32)
    tx, ty, tz = template_grid
    atlas = np.empty((*TEMPLATE_SHAPE, CLASSES), dtype=np.float32)
    for k in range(CLASSES):
        atlas[..., k] = np.clip(
            .22 + .09 * np.sin(.41 * tx + .13 * k) * np.cos(.37 * ty - .09 * k)
            + .04 * np.sin(.29 * tz + .17 * (tx + ty)), 0, 1
        )
    template_affine = AFFINE.copy()
    # Template voxel (1,1,1) occupies the same RAS point as subject voxel (0,0,0).
    template_affine[:3, 3] -= template_affine[:3, :3] @ np.ones(3)

    _write("images", image)
    _write("intensity_field", field)
    _write("psf_bias", bias)
    _write("initial_memberships", initial)
    _write("transformed_shapes", shapes)
    _write("closest_labels", closest)
    _write("folded_ribbon", folded.astype(np.uint8))
    _write("atlas_template_shapes", atlas, template_affine)


if __name__ == "__main__":
    main()
