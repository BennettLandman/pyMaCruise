from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "lesion_competition_pass"
SHAPE = (17, 17, 17)
CLASSES = 4


def create() -> np.ndarray:
    nx, ny, nz = SHAPE
    labels = np.zeros(SHAPE, dtype=np.int16)
    yy, zz = np.mgrid[:ny, :nz]
    # A curved, folded cortical ribbon with a short narrow branch. The two
    # folds vary independently so x/y/z scan order is observable.
    ridge = 8 + np.rint(1.5 * np.sin(yy * 0.55) + 1.1 * np.cos(zz * 0.47)).astype(int)
    for y in range(3, 14):
        for z in range(3, 14):
            center = int(ridge[y, z])
            labels[center - 2:center + 3, y, z] = 1
    labels[5:12, 8, 8] = 1
    labels[7:9, 7:9, 8:10] = 2
    labels[9, 9, 9] = 3

    mems = np.full((*SHAPE, CLASSES), 0.14, dtype=np.float32)
    # Spatially varying but deterministic probabilities ensure nonuniform
    # priorities and equal-priority heap ties in the competition frontier.
    xx, yy_full, zz_full = np.indices(SHAPE)
    for k in range(CLASSES):
        mems[..., k] += np.float32(0.015) * (
            (xx * (k + 1) + yy_full * (k + 2) + zz_full * (k + 3)) % 7
        ).astype(np.float32)
    for k in range(CLASSES):
        mems[..., k] += np.float32(labels == k) * np.float32(0.46)
    # In the broad background, non-background membership is zero. This makes
    # the selected source's six-face queue stop before a change could be
    # scheduled at the outer, non-topologized image border.
    background = labels == 0
    mems[background] = np.float32(0.0)
    # Competition's selected follow-on heap payload omits the target label;
    # background-class membership zero also makes its resulting 0->0 entries
    # drain instead of repeatedly scheduling self transitions.
    mems[..., 0] = np.float32(0.0)
    lesion = np.zeros(SHAPE, dtype=np.float32)
    lesion[7:9, 7:9, 8:10] = np.float32(0.23)
    lesion[9, 9, 9] = np.float32(0.09)
    blackhole = np.zeros(SHAPE, dtype=np.float32)
    blackhole[9, 9, 9] = np.float32(0.31)
    ordering = np.full(SHAPE, np.float32(-7.5), dtype=np.float32)
    packed = np.empty((*SHAPE, CLASSES + 4), dtype=np.float32)
    packed[..., 0] = labels
    packed[..., 1:1 + CLASSES] = mems
    packed[..., 1 + CLASSES] = lesion
    packed[..., 2 + CLASSES] = blackhole
    packed[..., 3 + CLASSES] = ordering
    return packed


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(create(), np.diag([0.8, 1.0, 1.2, 1.0]))
    image.header.set_xyzt_units("mm")
    image.header.set_slope_inter(1.0, 0.0)
    nib.save(image, ROOT / "competition_inputs.nii.gz")


if __name__ == "__main__":
    main()
