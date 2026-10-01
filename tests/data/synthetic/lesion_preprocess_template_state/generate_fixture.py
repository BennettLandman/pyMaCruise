from pathlib import Path

import nibabel as nib
import numpy as np


root = Path(__file__).parent
x, y, z = np.indices((5, 6, 7))
labels = (1 + ((x * 11 + y * 5 + z * 3 + x * y) % 13)).astype(np.uint8)
affine = np.diag((0.8, 1.0, 1.2, 1.0))
nib.save(nib.Nifti1Image(labels, affine), root / "transformed_template.nii.gz")
