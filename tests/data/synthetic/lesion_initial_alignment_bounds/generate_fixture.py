from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE_DIMS = (35, 37, 36)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])
shape_files: list[str] = []
shape_sha256: dict[str, str] = {}

for class_index in range(3):
    values = np.zeros(SHAPE_DIMS, dtype=np.float32)
    for x in range(SHAPE_DIMS[0]):
        for y in range(SHAPE_DIMS[1]):
            for z in range(SHAPE_DIMS[2]):
                if not (
                    3 + class_index <= x < 31 - class_index
                    and 4 <= y < 33 - class_index
                    and 2 + class_index <= z < 32
                ):
                    continue
                fold = (
                    3 * x
                    + 5 * y
                    + 7 * z
                    + 11 * x * y
                    + 13 * y * z
                    + class_index
                ) % 11 - 5
                value = np.float32(0.72 - 0.12 * class_index)
                value = np.float32(value + np.float32(0.0015 * x))
                value = np.float32(value + np.float32(0.0011 * y))
                value = np.float32(value + np.float32(0.0009 * z))
                values[x, y, z] = np.float32(
                    value + np.float32(fold * 0.0002)
                )
    filename = f"shape_{class_index}_bounded_ras.nii.gz"
    nib.save(nib.Nifti1Image(values, AFFINE), ROOT / filename)
    shape_files.append(filename)
    shape_sha256[filename] = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()

(ROOT / "inputs.json").write_text(
    json.dumps(
        {
            "source_alignment_fixture": "../lesion_initial_atlas_alignment",
            "shape_niftis": shape_files,
            "shape_sha256": shape_sha256,
            "shape_dimensions_xyz": list(SHAPE_DIMS),
            "axis_codes": list(nib.aff2axcodes(AFFINE)),
            "support_rule": "class k: x=[3+k,31-k), y=[4,33-k), z=[2+k,32)",
            "shape_formula": "selected folded positive values inside support; zero outside",
        },
        indent=2,
    )
    + "\n"
)
