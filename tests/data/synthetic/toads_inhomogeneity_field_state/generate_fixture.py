"""Generate the deterministic RAS NIfTI field used by the setter test."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
field = np.asarray(
    [
        [[1.0, 1.25], [0.75, 1.5]],
        [[0.8, 1.2], [1.3, 0.9]],
    ],
    dtype=np.float32,
)
affine = np.diag([0.8, 1.2, 1.7, 1.0])
nib.save(nib.Nifti1Image(field, affine), ROOT / "field.nii.gz")
