"""Regenerate the compact asymmetric folded TOADS iteration phantom."""
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).parent
shape = (13, 13, 11)
nx, ny, nz = shape
segmentation = np.zeros(shape, dtype=np.int16)
images = np.zeros((*shape, 1), dtype=np.float32)
memberships = np.zeros((*shape, 3), dtype=np.float32)
shapes = np.zeros((*shape, 3), dtype=np.float32)
for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            # A folded, oblique centerline with a one-sided notch keeps this
            # phantom asymmetric while providing cortex-like narrow bands.
            center = 5 + ((2 * y + z + (y // 3)) % 3) - 1
            dx = abs(x - center)
            if 2 <= y <= 10 and 1 <= z <= 9 and dx <= 1:
                label = 1
            elif 1 <= y <= 11 and 0 <= z <= 10 and dx <= 3 and not (y > 8 and z < 3):
                label = 2
            else:
                label = 0
            segmentation[x, y, z] = label
            image_value = np.float32(0.12 + 0.31 * label + 0.015 * x + 0.009 * y + 0.006 * z)
            images[x, y, z, 0] = image_value
            distances = np.asarray((0.0, abs(dx - 1.0), abs(dx - 3.0)), dtype=np.float32)
            scores = np.exp(-distances, dtype=np.float32) + np.asarray((0.07, 0.03, 0.02), dtype=np.float32)
            memberships[x, y, z] = scores / scores.sum(dtype=np.float32)
            shapes[x, y, z] = np.asarray(
                (0.15 + 0.01 * z, np.exp(-abs(dx - 1.0)), np.exp(-abs(dx - 3.0))),
                dtype=np.float32,
            )

# Array axes are declared as RAS-oriented x,y,z with positive R/A/S directions.
affine = np.diag((1.0, 1.25, 1.5, 1.0))
for name, data in (
    ("segmentation.nii.gz", segmentation),
    ("image_values.nii.gz", images),
    ("initial_memberships.nii.gz", memberships),
    ("transformed_shapes.nii.gz", shapes),
):
    nib.save(nib.Nifti1Image(data, affine), str(ROOT / name))
