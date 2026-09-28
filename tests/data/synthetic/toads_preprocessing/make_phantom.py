"""Generate the deterministic TOADS robust-normalization phantom."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (13, 11, 9)


def main() -> None:
    count = int(np.prod(SHAPE))
    linear = np.arange(count)
    x = linear % SHAPE[0]
    y = (linear // SHAPE[0]) % SHAPE[1]
    z = linear // (SHAPE[0] * SHAPE[1])
    mask = (
        (x >= 1)
        & (x < SHAPE[0] - 1)
        & (y >= 1)
        & (y < SHAPE[1] - 1)
        & (z >= 1)
        & (z < SHAPE[2] - 1)
    )
    values = (np.float32(0.5) + ((linear * 37) % 501).astype(np.float32) / 10).astype(
        np.float32
    )
    values[[351, 527, 703]] = np.asarray([651.0, 827.0, 1003.0], dtype=np.float32)
    values[~mask] = 0.0
    image = values.reshape(SHAPE, order="F")
    brain_mask = mask.astype(np.uint8).reshape(SHAPE, order="F")
    affine = np.eye(4, dtype=np.float64)
    nib.save(nib.Nifti1Image(image, affine), ROOT / "t1_n4.nii.gz")
    nib.save(nib.Nifti1Image(brain_mask, affine), ROOT / "brain_mask.nii.gz")
    (ROOT / "spec.json").write_text(
        json.dumps(
            {
                "shape": list(SHAPE),
                "affine": affine.tolist(),
                "foreground_voxels": int(mask.sum()),
                "masked_voxels": int((~mask).sum()),
                "recipe": "TOADS source-indexed ramp with three bright outliers and a rectangular brain mask",
                "java_reference": "tests/data/synthetic/toads_preprocessing/java_reference.json",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
