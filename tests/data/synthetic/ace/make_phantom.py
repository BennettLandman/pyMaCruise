"""Regenerate the committed, analytic RAS ACE input phantom."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
SHAPE = (24, 24, 24)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def main() -> None:
    grid = np.indices(SHAPE, dtype=np.float32)
    radius = np.sqrt(sum((axis - 11.5) ** 2 for axis in grid))
    white = np.where(radius <= 4.0, 230.0, 0.0).astype(np.float32)
    gray = np.where((radius > 5.0) & (radius <= 9.0), 220.0, 0.0).astype(np.float32)
    atlas_csf = np.where((radius > 5.0) & (radius <= 9.0), 254, 0).astype(np.int16)

    files = {
        "gray_matter.nii.gz": gray,
        "white_matter.nii.gz": white,
        "multi_atlas_csf.nii.gz": atlas_csf,
    }
    for name, data in files.items():
        nib.save(nib.Nifti1Image(data, AFFINE), ROOT / name)

    metadata = {
        "description": "Concentric, analytically defined ACE tissue phantom",
        "shape_voxels": SHAPE,
        "voxel_size_mm": [1.0, 1.0, 1.0],
        "orientation": "RAS",
        "center_voxel": [11.5, 11.5, 11.5],
        "white_matter": "230 for radius <= 4 mm; 0 elsewhere",
        "gray_matter": "220 for 5 mm < radius <= 9 mm; 0 elsewhere",
        "multi_atlas_csf": "254 for 5 mm < radius <= 9 mm; 0 elsewhere",
        "generation": "float32 tissue volumes; int16 multi-atlas CSF",
    }
    (ROOT / "spec.json").write_text(json.dumps(metadata, indent=2) + "\n")
    names = sorted([*files, "spec.json"])
    lines = []
    for name in names:
        digest = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        lines.append(f"{digest}  {name}")
    (ROOT / "SHA256SUMS").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
