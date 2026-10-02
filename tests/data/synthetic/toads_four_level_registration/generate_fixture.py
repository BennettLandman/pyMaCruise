"""Generate the small synthetic inputs for four-level TOADS registration."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent


def main() -> None:
    # Positive diagonal affine gives a canonical RAS grid. The subject is large
    # enough that startup padding and block averaging preserve 2+ voxels at 8x.
    subject_shape = (32, 32, 32)
    x, y, z = np.indices(subject_shape, dtype=np.float32)
    center = (np.asarray(subject_shape, dtype=np.float32) - 1) / 2
    radius = ((x - center[0]) / 12) ** 2 + ((y - center[1]) / 11) ** 2 + ((z - center[2]) / 10) ** 2
    brain_mask = (radius <= 1).astype(np.uint8)
    t1 = np.zeros(subject_shape, dtype=np.float32)
    inside = brain_mask > 0
    t1[inside] = (
        25.0
        + 0.63 * x[inside]
        + 0.29 * y[inside]
        + 0.17 * z[inside]
        + 0.11 * ((x[inside] + 2 * y[inside] + z[inside]) % 7)
    )

    # 16^3 is the minimum native atlas size that retains 2 samples per axis at
    # the coarsest 8x level; four asymmetric smooth priors make the optimizer's
    # transform and gradients nontrivial.
    atlas_shape = (16, 16, 16)
    ax, ay, az = np.indices(atlas_shape, dtype=np.float32)
    topology = np.ones(atlas_shape, dtype=np.uint8)
    folded = ((ax - 7.5) ** 2 + 1.2 * (ay - 7.5) ** 2 + 0.8 * (az - 7.5) ** 2) < 45 + 2 * np.sin(ax / 2)
    topology[folded] = 2 + ((ax[folded].astype(np.int32) + ay[folded].astype(np.int32) + az[folded].astype(np.int32)) % 3)
    centers = ((7.0, 7.5, 7.0), (5.5, 8.0, 7.5), (9.0, 6.5, 8.0), (7.5, 9.0, 5.5))
    atlas_shapes = np.empty(atlas_shape + (4,), dtype=np.float32)
    for channel, (cx, cy, cz) in enumerate(centers):
        distance = (ax - cx) ** 2 + 1.25 * (ay - cy) ** 2 + 0.8 * (az - cz) ** 2
        atlas_shapes[..., channel] = np.exp(-distance / np.float32(20.0))
    atlas_shapes[..., 0] = 1.0

    affine = np.diag((1.0, 1.0, 1.0, 1.0))
    ROOT.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(t1, affine), ROOT / "t1_ras.nii.gz")
    nib.save(nib.Nifti1Image(brain_mask, affine), ROOT / "brain_mask_ras.nii.gz")
    nib.save(nib.Nifti1Image(topology, affine), ROOT / "atlas_topology_ras.nii.gz")
    nib.save(nib.Nifti1Image(atlas_shapes, affine), ROOT / "atlas_shape_priors_ras.nii.gz")


if __name__ == "__main__":
    main()
