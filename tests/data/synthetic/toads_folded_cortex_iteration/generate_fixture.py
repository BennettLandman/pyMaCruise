"""Regenerate the larger deterministic folded-cortex TOADS phantom."""

from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).parent
shape = (27, 25, 23)
nx, ny, nz = shape
segmentation = np.zeros(shape, dtype=np.int16)
images = np.zeros((*shape, 1), dtype=np.float32)
memberships = np.zeros((*shape, 3), dtype=np.float32)
shapes = np.zeros((*shape, 3), dtype=np.float32)

for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            # An asymmetric, cortex-like closed band whose centerline folds
            # across both in-plane axes. Its width varies slowly with angle,
            # leaving narrow sulcal turns and a small one-sided notch.
            fold = 1.45 * np.sin((y + 0.4 * z) * 0.43) + 0.8 * np.sin((z - 3) * 0.57)
            cx = 13.0 + fold
            cy = 12.0 + 0.9 * np.sin((x - 5) * 0.24)
            cz = 11.0 + 0.65 * np.sin((y + 2) * 0.31)
            radius = np.sqrt(((x - cx) / 9.3) ** 2 + ((y - cy) / 8.1) ** 2 + ((z - cz) / 7.2) ** 2)
            local_fold = (np.sin(0.57 * y + 0.31 * z) + 1.0) * 0.5
            inner_edge = 0.77 + 0.035 * local_fold
            outer_edge = 1.00 + 0.045 * np.sin(0.29 * y - 0.37 * z)
            if inner_edge <= radius <= outer_edge and not (y > 18 and z < 6 and x > 17):
                label = 1  # thin cortical ribbon
            elif radius < outer_edge + 0.20:
                label = 2  # surrounding tissue class
            else:
                label = 0  # exterior/background class
            segmentation[x, y, z] = label

            distances = np.asarray(
                (max(0.0, outer_edge + 0.05 - radius), abs(radius - 0.89), abs(radius - 1.08)),
                dtype=np.float32,
            )
            logits = -distances * np.asarray((2.2, 2.8, 1.8), dtype=np.float32)
            logits += np.asarray((0.015 * x, 0.018 * y, 0.012 * z), dtype=np.float32)
            scores = np.exp(logits - logits.max(), dtype=np.float32)
            memberships[x, y, z] = scores / scores.sum(dtype=np.float32)
            shapes[x, y, z] = np.asarray(
                (max(0.0, outer_edge + 0.05 - radius), np.exp(-abs(radius - 0.89)), np.exp(-abs(radius - 1.08))),
                dtype=np.float32,
            )
            images[x, y, z, 0] = np.float32(0.22 + 0.19 * label + 0.012 * x + 0.008 * y + 0.005 * z + 0.025 * fold)

# Positive diagonal means RAS-oriented x/y/z array axes and 1 mm isotropic voxels.
affine = np.diag((1.0, 1.0, 1.0, 1.0))
for name, data in (
    ("segmentation.nii.gz", segmentation),
    ("image_values.nii.gz", images),
    ("initial_memberships.nii.gz", memberships),
    ("transformed_shapes.nii.gz", shapes),
):
    nib.save(nib.Nifti1Image(data, affine), str(ROOT / name))
