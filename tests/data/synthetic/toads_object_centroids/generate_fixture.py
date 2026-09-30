"""Regenerate the tiny asymmetric RAS NIfTI input for OBJECT-centroid tests."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
shape = (5, 4, 4)
image = np.zeros((*shape, 2), dtype=np.float32)
field = np.zeros_like(image)
memberships = np.zeros((*shape, 3), dtype=np.float32)
segmentation = np.zeros(shape, dtype=np.int16)
for x in range(shape[0]):
    for y in range(shape[1]):
        for z in range(shape[2]):
            segmentation[x, y, z] = (x + 2 * y + z) % 3
            for channel in range(2):
                image[x, y, z, channel] = np.float32(
                    np.float32(0.17)
                    + np.float32(
                        np.float32(((x * 7 + y * 5 + z * 3 + channel * 11) % 19)
                                   * np.float32(0.071))
                    )
                )
                field[x, y, z, channel] = np.float32(
                    np.float32(0.73)
                    + np.float32(
                        np.float32(((x * 3 + y * 2 + z * 5 + channel) % 7)
                                   * np.float32(0.047))
                    )
                )
            for class_index in range(3):
                memberships[x, y, z, class_index] = np.float32(
                    0.0
                    if class_index == 2
                    else np.float32(0.09)
                    + np.float32(
                        np.float32(((x * 2 + y * 3 + z * 5 + class_index * 7) % 11)
                                   * np.float32(0.063))
                    )
                )

affine = np.diag([1.0, 1.0, 1.0, 1.0])
nib.save(nib.Nifti1Image(np.concatenate((image, field, memberships), axis=3), affine), ROOT / "object_centroid_inputs.nii.gz")
nib.save(nib.Nifti1Image(segmentation, affine), ROOT / "segmentation.nii.gz")
