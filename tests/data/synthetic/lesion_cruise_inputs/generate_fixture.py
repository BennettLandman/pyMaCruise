"""Regenerate the small asymmetric NIfTI inputs for the Java excerpt harness."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUT = Path(__file__).parent
SHAPE = (5, 4, 3)
NAMES = (
    "Background", "SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM",
    "CerebralGM", "Cerebellum-GM", "Cerebellar-WM", "Brainstem", "Other",
    "Lesion", "BlackHole", "WhiteMatter",
)
AFFINE = np.asarray(
    [[0.8, 0.0, 0.0, 13.0], [0.0, 1.2, 0.0, -7.0], [0.0, 0.0, 2.5, 4.0], [0, 0, 0, 1]],
    dtype=np.float64,
)


def main() -> None:
    voxel = np.arange(np.prod(SHAPE), dtype=np.int32).reshape(SHAPE)
    memberships = np.empty((*SHAPE, len(NAMES)), dtype=np.float32)
    for cls in range(len(NAMES)):
        numerator = ((voxel + 1) * (cls + 3)) % 23
        memberships[..., cls] = numerator.astype(np.float32) / np.float32(23.0)
    lesion = ((voxel + 2) % 11).astype(np.float32) / np.float32(17.0)
    black_hole = ((voxel + 3) % 7).astype(np.float32) / np.float32(29.0)
    # The relatively prime multiplier walks through the label table while
    # retaining asymmetric x/y/z placements.
    segmentation = ((voxel * 7 + 5) % len(NAMES)).astype(np.int16)

    nib.save(nib.Nifti1Image(memberships, AFFINE), OUT / "memberships.nii.gz")
    nib.save(nib.Nifti1Image(lesion, AFFINE), OUT / "lesion_membership.nii.gz")
    nib.save(nib.Nifti1Image(black_hole, AFFINE), OUT / "black_hole_membership.nii.gz")
    nib.save(nib.Nifti1Image(segmentation, AFFINE), OUT / "segmentation.nii.gz")


if __name__ == "__main__":
    main()
