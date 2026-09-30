"""Generate the small true-RAS label phantoms for selected CSF distance."""

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import save_nifti_ras


HERE = Path(__file__).resolve().parent


def main() -> None:
    seeded = np.zeros((9, 10, 11), dtype=np.uint8)
    values = {
        (0, 2, 10): 1, (1, 3, 8): 2, (2, 4, 7): 3, (3, 5, 6): 4,
        (4, 6, 5): 5, (5, 7, 4): 6, (8, 9, 0): 7,
        (1, 8, 10): 8, (3, 0, 9): 9, (5, 2, 8): 10,
        (7, 4, 6): 11, (0, 9, 3): 12, (8, 0, 5): 13,
    }
    for point, label in values.items():
        seeded[point] = label

    near_misses = np.zeros((5, 6, 7), dtype=np.uint8)
    for point, label in {
        (0, 0, 0): 8, (2, 3, 5): 9, (4, 1, 6): 10,
        (3, 2, 1): 11, (1, 5, 3): 12, (4, 5, 0): 13,
    }.items():
        near_misses[point] = label

    affine = np.asarray(
        [[1.1, 0, 0, 17.0], [0, 1.7, 0, -23.0], [0, 0, 2.3, 8.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    save_nifti_ras(nib.Nifti1Image(seeded, affine), HERE / "seeded_atlas_indices.nii.gz", dtype=np.uint8)
    save_nifti_ras(nib.Nifti1Image(near_misses, affine), HERE / "near_miss_atlas_indices.nii.gz", dtype=np.uint8)


if __name__ == "__main__":
    main()
