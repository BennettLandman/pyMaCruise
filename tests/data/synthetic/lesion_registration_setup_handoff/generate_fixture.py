from pathlib import Path

import nibabel as nib
import numpy as np


root = Path(__file__).parent
memberships = np.zeros((5, 6, 7, 3), dtype=np.float32)
affine = np.diag((0.8, 1.0, 1.2, 1.0))
nib.save(nib.Nifti1Image(memberships, affine), root / "initial_memberships.nii.gz")
