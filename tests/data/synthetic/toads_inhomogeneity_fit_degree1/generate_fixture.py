"""Generate a masked top-k fit with a known linear multiplicative field."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
AFFINE = np.diag([.85, 1.15, 1.55, 1.0])
DIMENSIONS = (9, 9, 9)


def main() -> None:
    x, y, z = np.indices(DIMENSIONS, dtype=np.float32)
    xn, yn, zn = x / np.float32(8), y / np.float32(8), z / np.float32(8)
    known_field = np.float32(1.05) + np.float32(.2) * xn - np.float32(.11) * yn + np.float32(.07) * zn
    image = np.float32(.72) / known_field
    memberships = np.empty(DIMENSIONS + (3,), dtype=np.float32)
    memberships[..., 0] = np.float32(.31) + np.float32(.02) * ((x + y) % 4)
    memberships[..., 1] = np.float32(.42) + np.float32(.025) * ((y + 2 * z) % 5)
    memberships[..., 2] = np.float32(.28) + np.float32(.03) * ((2 * x + z) % 4)
    best = np.empty(DIMENSIONS + (3,), dtype=np.uint8)
    best[..., 0] = 0  # source mask class is excluded from the polynomial fit
    best[..., 1] = 1
    best[..., 2] = 2
    segmentation = ((x + 2 * y + z) % 7 != 0).astype(np.uint8)
    for name, data in (("image", image), ("memberships", memberships), ("best", best), ("segmentation", segmentation)):
        nii = nib.Nifti1Image(data, AFFINE)
        nii.set_qform(AFFINE, code=1)
        nii.set_sform(AFFINE, code=1)
        nib.save(nii, ROOT / f"{name}.nii.gz")
    centers = np.full(11, np.float32(.72), dtype=np.float32)
    centers[0] = np.float32(.15)
    np.savetxt(ROOT / "centroids.txt", centers, fmt="%.9g")


if __name__ == "__main__":
    main()
