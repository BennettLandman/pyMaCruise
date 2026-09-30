from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


root = Path(__file__).parent
reference = json.loads((root / "java_reference.json").read_text())
shape = tuple(reference["shape"])
affine = np.eye(4, dtype=np.float64)
affine[:3, :3] = np.diag([1.0, 1.0, 1.0])

arrays = {
    "classification": np.asarray(reference["classification"], dtype=np.int8).reshape(shape),
    "memberships": np.asarray(reference["memberships"], dtype=np.float32).reshape(
        (*shape, len(reference["template_labels"]))
    ),
    "lesion_membership": np.asarray(reference["lesion_membership"], dtype=np.float32).reshape(shape),
}
for name, array in arrays.items():
    image = nib.Nifti1Image(array, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, root / f"{name}.nii.gz")
