"""Generate a deterministic 41³ radially folded synthetic cortex phantom."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (41, 41, 41)
CENTER = np.array([20.0, 20.0, 20.0], dtype=np.float32)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def phantom_arrays() -> dict[str, np.ndarray]:
    x, y, z = np.indices(SHAPE, dtype=np.float32)
    dx, dy, dz = x - CENTER[0], y - CENTER[1], z - CENTER[2]
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    theta = np.arctan2(dy, dx)
    phi = np.arctan2(np.sqrt(dx * dx + dy * dy), dz)
    # Angular lobes model repeated cortical folds while keeping nested,
    # star-shaped tissue surfaces for deterministic surface/narrow-band tests.
    fold = 2.0 * np.sin(phi) ** 2 * np.cos(6.0 * theta)
    fold += 1.0 * np.sin(phi) ** 2 * np.cos(4.0 * phi)
    outer_radius = 13.0 + fold
    inner_radius = outer_radius - 4.5
    inside = radius <= outer_radius
    white = radius <= inner_radius
    gray = inside & ~white
    mid_surface = np.abs(radius - (outer_radius + inner_radius) * 0.5) <= 0.55
    return {
        "initial_membership": np.where(inside, 255, 0).astype(np.uint8),
        "white_matter": np.where(white, 230, 0).astype(np.float32),
        "gray_matter": np.where(gray, 220, 0).astype(np.float32),
        "ace_gray_matter": np.where(gray, 220, 0).astype(np.float32),
        "multi_atlas_csf": np.where(gray, 254, 0).astype(np.int16),
        "ace_skeleton": np.where(mid_surface, 254, 0).astype(np.uint8),
    }


def main() -> None:
    arrays = phantom_arrays()
    for name, values in arrays.items():
        image = nib.Nifti1Image(values, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")

    description = {
        "description": "Synthetic cortex-like surface with deterministic angular folds; no subject data.",
        "shape_voxels": list(SHAPE),
        "center_voxel": CENTER.tolist(),
        "voxel_size_mm": [1.0, 1.0, 1.0],
        "orientation": "RAS",
        "outer_radius_mm": "13 + 2*sin(phi)^2*cos(6*theta) + sin(phi)^2*cos(4*phi)",
        "cortical_thickness_mm": 4.5,
        "generation": "deterministic analytic radial formula; NumPy float32 coordinates",
        "roles": {
            "initial_membership": "255 inside the folded outer surface",
            "white_matter": "230 inside the folded inner surface",
            "gray_matter": "220 between inner and outer surfaces",
            "multi_atlas_csf": "254 in the cortical shell",
            "ace_skeleton": "254 within 0.55 mm of the shell mid-surface",
        },
    }
    (ROOT / "spec.json").write_text(json.dumps(description, indent=2) + "\n")
    payloads = sorted(ROOT.glob("*.nii.gz")) + [ROOT / "spec.json"]
    (ROOT / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
            for path in payloads
        )
    )


if __name__ == "__main__":
    main()

