from pathlib import Path

import nibabel as nib
import numpy as np


root = Path(__file__).parent
x, y, z = np.indices((9, 10, 11))
labels = (1 + ((x * 7 + y * 3 + z * 5 + x * y + 2 * z * y) % 17)).astype(np.uint8)
labels[0, 0, 0] = 42  # differs from the selected atlas metadata label[0] = 1
affine = np.diag((0.9, 1.0, 1.1, 1.0))
nib.save(nib.Nifti1Image(labels, affine), root / "template_labels.nii.gz")
