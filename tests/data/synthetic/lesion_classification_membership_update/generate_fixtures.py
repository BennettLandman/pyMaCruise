from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


root = Path(__file__).parent
reference = json.loads((root / "java_reference.json").read_text())
affine = np.eye(4, dtype=np.float64)
affine[:3, :3] = np.diag([1.0, 1.0, 1.0])

for case in reference["scenarios"]:
    shape = tuple(case["shape"])
    arrays = {
        "previous": np.asarray(case["previous"], dtype=np.int8).reshape(shape),
        "memberships": np.asarray(case["memberships"], dtype=np.float32).reshape(
            (*shape, len(case["template_labels"]))
        ),
        "lesion_membership": np.asarray(case["lesion_membership"], dtype=np.float32).reshape(shape),
    }
    for name, array in arrays.items():
        image = nib.Nifti1Image(array, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, root / f"{case['name']}_{name}.nii.gz")
