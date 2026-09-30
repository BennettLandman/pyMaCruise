from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUT = Path(__file__).parent
NAMES = (
    "SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM",
    "CerebralGM", "Cerebellum-GM", "CerebellarGM", "CerebellarGM)",
    "WhiteMatter", "Background", "Other", "other",
)
SHAPE = (5, 2, 2)
AFFINE = np.asarray(
    [[1.1, 0.0, 0.0, 12.0], [0.0, 1.2, 0.0, -23.0], [0.0, 0.0, 1.3, 34.0], [0, 0, 0, 1]],
    dtype=np.float64,
)


def inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x, y, z = np.indices(SHAPE, dtype=np.int32)
    labels = ((x + 2 * y + 3 * z) % len(NAMES)).astype(np.int16)
    lesion = (0.013 * (1 + x + 3 * y + 7 * z)).astype(np.float32)
    black = (0.017 * (2 + 2 * x + y + 5 * z)).astype(np.float32)
    memberships = np.empty((*SHAPE, len(NAMES)), dtype=np.float32)
    for k in range(len(NAMES)):
        memberships[..., k] = (0.001 * (1 + x + 5 * y + 11 * z + 13 * k)).astype(np.float32)
    return memberships, labels, lesion, black


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, array in zip(
        ("memberships.nii.gz", "segmentation.nii.gz", "lesion_membership.nii.gz", "black_hole_membership.nii.gz"),
        inputs(),
    ):
        nib.save(nib.Nifti1Image(array, AFFINE), OUT / name)


if __name__ == "__main__":
    main()
