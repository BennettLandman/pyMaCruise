"""Generate a compact known-transform RAS registration fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import ants


OUTPUT = Path(__file__).parent
SHAPE = (17, 17, 17)
AFFINE = np.eye(4)
MOVING_CENTER = np.array([5.0, 8.0, 7.0])
SHIFT_RAS_MM = np.array([2.0, 0.0, 0.0])


def phantom(center: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grid = np.indices(SHAPE, dtype=np.float32)
    squared_distance = sum((grid[axis] - center[axis]) ** 2 for axis in range(3))
    gaussian = np.exp(-squared_distance / 5.0).astype(np.float32)
    labels = np.zeros(SHAPE, dtype=np.uint8)
    labels[squared_distance <= 2.5**2] = 7
    return gaussian, labels


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    moving_image, moving_labels = phantom(MOVING_CENTER)
    moving_flair = np.sqrt(moving_image).astype(np.float32)
    fixed_image, fixed_labels = phantom(MOVING_CENTER + SHIFT_RAS_MM)
    mask = np.zeros(SHAPE, dtype=np.uint8)
    mask[2:-2, 2:-2, 2:-2] = 1
    for name, values in {
        "fixed_image": fixed_image,
        "moving_image": moving_image,
        "moving_flair": moving_flair,
        "fixed_labels": fixed_labels,
        "moving_labels": moving_labels,
        "fixed_mask": mask,
    }.items():
        nib.save(nib.Nifti1Image(values, AFFINE), str(OUTPUT / f"{name}.nii.gz"))
    # ANTs image resampling uses a fixed-grid-to-moving-image sampling map.
    # Therefore the file applied to map the moving sphere into the fixed image
    # has the opposite direction from moving-to-fixed point coordinates.
    forward = ants.create_ants_transform(
        transform_type="AffineTransform",
        precision="double",
        dimension=3,
        matrix=np.eye(3),
        translation=(2.0, 0.0, 0.0),
        center=(0.0, 0.0, 0.0),
    )
    inverse = ants.create_ants_transform(
        transform_type="AffineTransform",
        precision="double",
        dimension=3,
        matrix=np.eye(3),
        translation=(-2.0, 0.0, 0.0),
        center=(0.0, 0.0, 0.0),
    )
    ants.write_transform(forward, str(OUTPUT / "known_moving_to_fixed.mat"))
    ants.write_transform(inverse, str(OUTPUT / "known_fixed_to_moving.mat"))
    spec = {
        "purpose": "Synthetic ANTs transform direction/interpolation validation",
        "shape_voxels": list(SHAPE),
        "spacing_mm": [1, 1, 1],
        "orientation": "RAS",
        "moving_landmark_ras_mm": MOVING_CENTER.tolist(),
        "fixed_landmark_ras_mm": (MOVING_CENTER + SHIFT_RAS_MM).tolist(),
        "moving_to_fixed_ras_translation_mm": SHIFT_RAS_MM.tolist(),
        "ants_forward_image_resampling_lps_translation_mm": [2, 0, 0],
        "moving_to_fixed_point_translation_ras_mm": [2, 0, 0],
        "labels": {"background": 0, "sphere": 7},
        "seed": None,
    }
    (OUTPUT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    names = [
        "fixed_image.nii.gz",
        "moving_image.nii.gz",
        "moving_flair.nii.gz",
        "fixed_labels.nii.gz",
        "moving_labels.nii.gz",
        "fixed_mask.nii.gz",
        "known_moving_to_fixed.mat",
        "known_fixed_to_moving.mat",
        "spec.json",
    ]
    entries = [f"{hashlib.sha256((OUTPUT / name).read_bytes()).hexdigest()}  {name}" for name in names]
    (OUTPUT / "SHA256SUMS").write_text("\n".join(entries) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
