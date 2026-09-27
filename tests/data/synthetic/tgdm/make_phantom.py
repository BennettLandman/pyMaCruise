"""Generate a compact concentric RAS NIfTI phantom for TGDM tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (25, 25, 25)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def main() -> None:
    x, y, z = np.indices(SHAPE, dtype=np.float32)
    radius = np.sqrt((x - 12) ** 2 + (y - 12) ** 2 + (z - 12) ** 2)
    arrays = {
        "initial_membership": np.where(radius <= 8, 255, 0).astype(np.uint8),
        "white_matter": np.where(radius <= 5, 255, 0).astype(np.uint8),
        "gray_matter": np.where((radius > 5) & (radius <= 8), 220, 0).astype(np.uint8),
        "ace_gray_matter": np.where((radius > 5) & (radius <= 8), 220, 0).astype(np.uint8),
        "ace_skeleton": np.where(
            (np.abs(radius - 6.5) < 0.5) & (z == 12), 254, 0
        ).astype(np.uint8),
    }
    for name, values in arrays.items():
        image = nib.Nifti1Image(values, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, str(ROOT / f"{name}.nii.gz"))

    (ROOT / "spec.json").write_text(
        json.dumps(
            {
                "description": "Synthetic nested spherical tissue phantom; no subject data.",
                "shape": SHAPE,
                "affine": AFFINE.tolist(),
                "orientation": "RAS",
                "voxel_size_mm": [1, 1, 1],
                "center_voxel": [12, 12, 12],
                "radii_voxels": {"white_matter": 5, "gray_matter_outer": 8},
                "membership_values": [0, 255],
                "intensity_values": [0, 220, 255],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    records = []
    for path in sorted(ROOT.glob("*.nii.gz")):
        records.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (ROOT / "SHA256SUMS").write_text("\n".join(records) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
