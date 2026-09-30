"""Generate the folded RAS P18czl synthetic fixed-input fixture."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_lesion_atlas_type_membership"
SUBJECT = (9, 8, 7)
TEMPLATE = (11, 10, 9)
CLASSES, CHANNELS, LESION_SLOT = 13, 2, 10
SUBJECT_AFFINE = np.asarray(
    [[.9, 0, 0, -31.25], [0, 1.2, 0, -42.5], [0, 0, 1.45, 18.75], [0, 0, 0, 1]],
    dtype=np.float32,
)


def _save(name: str, array: np.ndarray, affine: np.ndarray = SUBJECT_AFFINE) -> None:
    image = nib.Nifti1Image(np.asarray(array), affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, str(ROOT / f"{name}.nii.gz"))


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    x, y, z = np.indices(SUBJECT, dtype=np.float32)
    fold = (x - 4) ** 2 / 9 + (y - 3) ** 2 / 5 + (z - 3) ** 2 / 4 < 1
    fold ^= ((x + 2 * y + z) % 5 == 0)
    images = np.empty((*SUBJECT, CHANNELS), dtype=np.float32)
    images[..., 0] = .16 + .029 * x + .018 * y + .014 * z + .063 * fold
    images[..., 1] = .12 + .015 * x + .024 * y + .019 * z + .087 * (fold & ((x + z) % 2 == 0))
    images[0, 2, 3, 0] = np.float32(0.0)
    images[-1, 3, 2, 1] = np.float32(-0.0)
    field = np.empty_like(images)
    bias = np.empty_like(images)
    for c in range(CHANNELS):
        field[..., c] = .92 + .005 * x + .004 * y + .003 * z + .012 * c
        bias[..., c] = .0018 * (x - 2) - .0011 * (y + 1) + .0013 * z + .0006 * c

    memberships = np.empty((*SUBJECT, CLASSES), dtype=np.float32)
    shapes = np.empty_like(memberships)
    for k in range(CLASSES):
        memberships[..., k] = .031 + .0025 * ((x + 2 * y + z + 2 * k) % 7)
        shapes[..., k] = .2 + .011 * ((2 * x + y + 3 * z + k) % 9)

    closest = np.empty((*SUBJECT, 4), dtype=np.int16)
    segmentation = np.empty(SUBJECT, dtype=np.int16)
    lesion_membership = np.empty(SUBJECT, dtype=np.float32)
    for ix in range(SUBJECT[0]):
        for iy in range(SUBJECT[1]):
            for iz in range(SUBJECT[2]):
                closest[ix, iy, iz] = (
                    (ix + 2 * iy + iz) % CLASSES,
                    (ix + iy + 3 * iz + 4) % CLASSES,
                    (2 * ix + iy + iz + 1) % CLASSES,
                    (ix + 3 * iy + iz + 7) % CLASSES,
                )
                if (ix + iy + iz) % 3 == 0:
                    segmentation[ix, iy, iz] = LESION_SLOT
                else:
                    segmentation[ix, iy, iz] = (ix + 2 * iy + 3 * iz) % CLASSES
                lesion_membership[ix, iy, iz] = .11 + .025 * ((ix + 2 * iy + iz) % 6)

    tx, ty, tz = np.indices(TEMPLATE, dtype=np.float32)
    template = np.empty((*TEMPLATE, CLASSES), dtype=np.float32)
    for k in range(CLASSES):
        template[..., k] = np.clip(
            .23 + .08 * np.sin(.39 * tx + .11 * k) * np.cos(.34 * ty - .07 * k)
            + .035 * np.sin(.27 * tz + .15 * (tx + ty)), 0, 1
        )
    template_affine = SUBJECT_AFFINE.copy()
    template_affine[:3, 3] -= template_affine[:3, :3] @ np.ones(3, dtype=np.float32)

    _save("images", images)
    _save("intensity_field", field)
    _save("psf_bias", bias)
    _save("initial_memberships", memberships)
    _save("transformed_shapes", shapes)
    _save("closest_labels", closest)
    _save("segmentation", segmentation)
    _save("initial_lesion_membership", lesion_membership)
    _save("folded_ribbon", fold.astype(np.uint8))
    _save("atlas_template_shapes", template, template_affine)


if __name__ == "__main__":
    main()
