"""Write deterministic true-RAS inputs for the selected 11-class chain."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> None:
    root = Path(__file__).parent
    size, classes = 14, 11
    x, y, z = np.indices((size, size, size), dtype=np.int16)
    labels = ((x * 3 + y * 5 + z * 7 + (x * y) % classes) % classes).astype(np.int16)
    skeleton = np.full(labels.shape, -1, dtype=np.int16)
    for k in range(classes):
        i, j, l = 1 + k, 1 + (k * 3) % 12, 1 + (k * 5) % 12
        labels[i, j, l] = k
        skeleton[i, j, l] = k

    memberships = np.empty((*labels.shape, classes), dtype=np.float32)
    for k in range(classes):
        memberships[..., k] = np.where(labels == k, np.float32(0.91), np.float32(0.009))
    ordering = -((x * 13 + y * 7 + z * 3) % 19).astype(np.float32) / np.float32(10.0)
    shapes = np.zeros((*labels.shape, classes), dtype=np.float32)
    affine = np.diag([0.8, 1.2, 1.7, 1.0]).astype(np.float64)
    affine[:3, 3] = [-12.0, 3.5, 8.0]
    for name, values in (
        ("initial_segmentation.nii.gz", labels),
        ("initial_skeleton.nii.gz", skeleton),
        ("exact_memberships.nii.gz", memberships),
        ("growing_ordering.nii.gz", ordering),
        ("transformed_shapes.nii.gz", shapes),
    ):
        nib.save(nib.Nifti1Image(values, affine), str(root / name))


if __name__ == "__main__":
    main()
