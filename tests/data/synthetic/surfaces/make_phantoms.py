"""Generate deterministic analytic RAS level-set fixtures for P30."""

from __future__ import annotations

import json
import hashlib
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
    thickness_output = OUTPUT / "thickness"
    thickness_output.mkdir(parents=True, exist_ok=True)
    for name, radius in {"inner": 5.5, "central": 6.5, "outer": 7.5}.items():
        # Enlarge the synthetic grid so all three nested spheres are enclosed.
        thickness_shape = (41, 41, 41)
        affine = np.diag([*SPACING_MM, 1.0])
        center_world = affine[:3, :3] @ ((np.asarray(thickness_shape) - 1) / 2)
        indices = np.indices(thickness_shape, dtype=np.float32)
        world = np.einsum("ij,j...->i...", affine[:3, :3], indices)
        values = np.linalg.norm(world - center_world.reshape(3, 1, 1, 1), axis=0) - radius
        thickness_image = nib.Nifti1Image(values.astype(np.float32), affine)
        nib.save(thickness_image, str(thickness_output / f"{name}.nii.gz"))
    (OUTPUT / "spec.json").write_text(
        json.dumps(
            {
                "purpose": "P30 synthetic level-set surface extraction validation",
                "shape_voxels": list(SHAPE),
                "spacing_mm": list(SPACING_MM),
                "thickness_shape_voxels": [41, 41, 41],
                "orientation": "RAS",
                "center_voxel": [(value - 1) / 2 for value in SHAPE],
                "signed_field": "min(radii) * (sqrt(sum(((world-center)/radii)^2)) - 1)",
                "radii_mm": {"sphere": [5.5, 5.5, 5.5], "ellipsoid": [5.5, 4.5, 6.0]},
                "thickness_sphere_radii_mm": {"inner": 5.5, "central": 6.5, "outer": 7.5},
                "seed": None,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    for directory, names in (
        (OUTPUT, ("sphere.nii.gz", "ellipsoid.nii.gz", "spec.json")),
        (
            thickness_output,
            ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"),
        ),
    ):
        entries = []
        for name in names:
            path = directory / name
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            entries.append(f"{digest}  {name}")
        (directory / "SHA256SUMS").write_text("\n".join(entries) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
