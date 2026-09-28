"""Regenerate the deterministic two-channel LesionTOADS initialization fixture."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (9, 8, 7)
LESION_PRIORS = np.asarray([0.4, 1.2], dtype=np.float32)
BLACK_HOLE_PRIORS = np.asarray([0.1, 0.6], dtype=np.float32)


def make_phantom() -> np.ndarray:
    count = int(np.prod(SHAPE))
    linear = np.arange(count, dtype=np.int32)
    images = np.empty((2, count), dtype=np.float32)
    images[0] = np.float32(4.25) + ((linear * 19 + 3) % 137).astype(np.float32) / np.float32(16.0)
    images[1] = np.float32(11.5) + ((linear * 29 + 5) % 211).astype(np.float32) / np.float32(12.0)
    images[:, ::61] = np.stack(
        (
            np.asarray([-5.0 - i for i in linear[::61]], dtype=np.float32),
            np.asarray([-9.0 - i for i in linear[::61]], dtype=np.float32),
        )
    )
    images[0, 37] = np.float32(600.0)
    images[0, 414] = np.float32(1200.0)
    images[1, 73] = np.float32(900.0)
    images[1, 422] = np.float32(1600.0)
    return images.reshape((2, *SHAPE), order="F")


def main() -> None:
    images = make_phantom()
    affine = np.eye(4, dtype=np.float64)
    nib.save(nib.Nifti1Image(images[0], affine), ROOT / "t1_mprage.nii.gz")
    nib.save(nib.Nifti1Image(images[1], affine), ROOT / "flair.nii.gz")
    (ROOT / "spec.json").write_text(
        json.dumps(
            {
                "shape_xyz": list(SHAPE),
                "channels": ["T1_MPRAGE", "FLAIR"],
                "affine_ras_mm": affine.tolist(),
                "robust_fraction": 0.01,
                "robust_scales": 2,
                "lesion_priors": LESION_PRIORS.tolist(),
                "black_hole_priors": BLACK_HOLE_PRIORS.tolist(),
                "recipe": (
                    "Two source-indexed modular ramps with periodic negative values, "
                    "and four high outliers; negatives are clamped to zero before "
                    "source robust bounds are computed."
                ),
                "java_reference": "java_reference.json",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
