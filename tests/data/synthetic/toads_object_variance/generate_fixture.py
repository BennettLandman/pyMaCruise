from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def make_fixture() -> None:
    x, y, z = np.indices((4, 3, 3), dtype=np.int32)
    images = np.empty((4, 3, 3, 2), dtype=np.float32)
    field = np.empty_like(images)
    memberships = np.empty((4, 3, 3, 4), dtype=np.float32)
    for c in range(2):
        image_step = ((7 * x + 5 * y + 3 * z + 11 * c) % 23).astype(np.float32)
        field_step = ((3 * x + 2 * y + 5 * z + c) % 9).astype(np.float32)
        images[..., c] = np.float32(0.13) + image_step * np.float32(0.017)
        field[..., c] = np.float32(0.71) + field_step * np.float32(0.031)
    for k in range(3):
        if k == 1:
            memberships[..., k] = np.float32(0.0)
        else:
            mem_step = ((3 * x + y + 2 * z + k) % 5).astype(np.float32)
            memberships[..., k] = mem_step / np.float32(4.0)
    segmentation = ((x + y + z) % 2).astype(np.int16)
    segmentation[0, 0, 0] = 2
    segmentation[-1, -1, -1] = 2
    for name, data in (
        ("images.nii.gz", images),
        ("intensity_field.nii.gz", field),
        ("memberships.nii.gz", memberships),
        ("segmentation.nii.gz", segmentation),
    ):
        nib.save(nib.Nifti1Image(data, AFFINE), str(ROOT / name))


if __name__ == "__main__":
    make_fixture()
