"""Regenerate analytic scalar-field inputs for the GVF port."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
SHAPE = (16, 12, 10)


def main() -> None:
    ramp = np.broadcast_to(
        np.arange(SHAPE[0], dtype=np.float32)[:, None, None], SHAPE
    ).copy()
    nib.save(nib.Nifti1Image(ramp, np.eye(4)), ROOT / "gradient_ramp.nii.gz")
    spec = {
        "description": "Linear scalar ramp with a known constant GVF gradient",
        "file": "gradient_ramp.nii.gz",
        "shape_voxels": SHAPE,
        "voxel_size_mm": [1.0, 1.0, 1.0],
        "orientation": "RAS",
        "values": "voxel index 0 along axis 0; constant along axes 1 and 2",
        "expected_normalized_gradient": [1.0 / (SHAPE[0] - 1), 0.0, 0.0],
    }
    (ROOT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    names = ["gradient_ramp.nii.gz", "spec.json"]
    (ROOT / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256((ROOT / name).read_bytes()).hexdigest()}  {name}\n"
            for name in names
        )
    )


if __name__ == "__main__":
    main()
