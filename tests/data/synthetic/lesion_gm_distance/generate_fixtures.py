"""Create the small true-RAS label phantom used by the GM-distance excerpt."""

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import save_nifti_ras


HERE = Path(__file__).resolve().parent


def main() -> None:
    seeded = np.zeros((9, 10, 11), dtype=np.uint8)
    seeds = {
        (0, 2, 10): 1, (1, 3, 8): 2, (2, 4, 7): 3, (3, 5, 6): 4,
        (4, 6, 5): 5, (5, 7, 4): 6, (6, 8, 3): 7, (7, 1, 2): 8,
        (8, 9, 0): 9, (2, 1, 1): 10, (4, 2, 9): 11,
        (6, 3, 7): 12, (8, 0, 5): 13,
    }
    controls = {
        (1, 8, 10): 14, (3, 0, 9): 15, (5, 2, 8): 16,
        (7, 4, 6): 17, (0, 9, 3): 18,
    }
    for point, label in (seeds | controls).items():
        seeded[point] = label

    aliases = np.zeros((5, 6, 7), dtype=np.uint8)
    for point, label in {
        (0, 0, 0): 14,
        (2, 3, 5): 15,
        (4, 1, 6): 16,
        (3, 2, 1): 17,
        (1, 5, 3): 18,
    }.items():
        aliases[point] = label

    # Positive RAS axes, anisotropic spacing, and nonzero origin.
    affine = np.asarray(
        [[1.1, 0, 0, 17.0], [0, 1.7, 0, -23.0], [0, 0, 2.3, 8.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    save_nifti_ras(nib.Nifti1Image(seeded, affine), HERE / "seeded_atlas_indices.nii.gz", dtype=np.uint8)
    save_nifti_ras(nib.Nifti1Image(aliases, affine), HERE / "near_miss_atlas_indices.nii.gz", dtype=np.uint8)


if __name__ == "__main__":
    main()
