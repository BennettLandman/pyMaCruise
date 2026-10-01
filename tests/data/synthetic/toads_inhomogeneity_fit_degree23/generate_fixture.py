"""Generate the folded anisotropic RAS phantom for degree-2/3 field fitting."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
shape = (18, 18, 18)
x, y, z = np.indices(shape, dtype=np.float32)
fold = (
    np.float32(7.2)
    + np.float32(1.8) * np.sin(np.float32(2.0 * np.pi / 9.0) * y)
    + np.float32(1.2) * np.cos(np.float32(2.0 * np.pi / 10.0) * z)
    + np.float32(0.35) * np.sin(np.float32(2.0 * np.pi / 6.0) * (y + z))
)
mask = (np.abs(x - fold) < np.float32(5.5)) & (y >= 1) & (y <= 16) & (z >= 1) & (z <= 16)
segmentation = np.where(mask, 1 + ((x + 2 * y + z).astype(np.int32) % 10), 0).astype(np.uint8)
image = np.asarray(
    np.float32(0.24)
    + np.float32(0.004) * x
    + np.float32(0.003) * y
    + np.float32(0.002) * z
    + np.float32(0.025) * np.sin(np.float32(2.0 * np.pi / 13.0) * (y + np.float32(0.4) * z))
    + np.float32(0.01) * np.cos(np.float32(2.0 * np.pi / 11.0) * (x + z)),
    dtype=np.float32,
)
memberships = np.zeros(shape + (3,), dtype=np.float32)
memberships[..., 0] = np.float32(0.1)
memberships[..., 1] = np.asarray(
    np.float32(0.3) + np.float32(0.4) * (np.float32(0.5) + np.float32(0.5) * np.sin((x + y) * np.float32(0.3))),
    dtype=np.float32,
)
memberships[..., 2] = np.asarray(
    np.float32(0.2) + np.float32(0.5) * (np.float32(0.5) + np.float32(0.5) * np.cos((y + z) * np.float32(0.27))),
    dtype=np.float32,
)
best = np.zeros(shape + (3,), dtype=np.uint8)
best[..., 1] = (1 + (x + y + z).astype(np.int32) % 10).astype(np.uint8)
best[..., 2] = (1 + (3 * x + 2 * y + z + 1).astype(np.int32) % 10).astype(np.uint8)
centroids = np.linspace(0.20, 0.56, 11, dtype=np.float32)
affine = np.diag([0.85, 1.15, 1.55, 1.0])

for name, array in (
    ("image", image),
    ("memberships", memberships),
    ("best", best),
    ("segmentation", segmentation),
):
    nib.save(nib.Nifti1Image(array, affine), ROOT / f"{name}.nii.gz")
np.savetxt(ROOT / "centroids.txt", centroids, fmt="%.9g")
