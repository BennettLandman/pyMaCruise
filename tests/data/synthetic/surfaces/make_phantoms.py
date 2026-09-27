"""Generate deterministic analytic RAS level-set fixtures for P30."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parent
SHAPE = (31, 31, 31)
SPACING_MM = (0.5, 0.5, 0.5)


def make_levelset(radii_mm: tuple[float, float, float]) -> nib.Nifti1Image:
    affine = np.diag([*SPACING_MM, 1.0])
    center_index = (np.asarray(SHAPE) - 1) / 2
    center_world = affine[:3, :3] @ center_index
    indices = np.indices(SHAPE, dtype=np.float32)
    world = np.einsum("ij,j...->i...", affine[:3, :3], indices)
    radii = np.asarray(radii_mm, dtype=np.float32).reshape(3, 1, 1, 1)
    normalized = (world - center_world.reshape(3, 1, 1, 1)) / radii
    values = (np.sqrt(np.sum(normalized**2, axis=0)) - 1.0) * min(radii_mm)
    return nib.Nifti1Image(values.astype(np.float32), affine)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, radii in {
        "sphere": (5.5, 5.5, 5.5),
        "ellipsoid": (5.5, 4.5, 6.0),
    }.items():
        nib.save(make_levelset(radii), str(OUTPUT / f"{name}.nii.gz"))
    (OUTPUT / "spec.json").write_text(
        json.dumps(
            {
                "purpose": "P30 synthetic level-set surface extraction validation",
                "shape_voxels": list(SHAPE),
                "spacing_mm": list(SPACING_MM),
                "orientation": "RAS",
                "center_voxel": [(value - 1) / 2 for value in SHAPE],
                "signed_field": "min(radii) * (sqrt(sum(((world-center)/radii)^2)) - 1)",
                "radii_mm": {"sphere": [5.5, 5.5, 5.5], "ellipsoid": [5.5, 4.5, 6.0]},
                "seed": None,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
