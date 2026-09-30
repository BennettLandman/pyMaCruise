"""Create the small, asymmetric, RAS-oriented NIfTI inputs for this slice."""

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import save_nifti_ras


HERE = Path(__file__).resolve().parent


def _case(shape: tuple[int, int, int], points: dict[tuple[int, int, int], int]) -> np.ndarray:
    labels = np.zeros(shape, dtype=np.uint8)
    for point, label in points.items():
        labels[point] = label
    return labels


def main() -> None:
    seeded = _case(
        (7, 8, 9),
        {
            (0, 2, 7): 1,  # border Brainstem
            (4, 1, 3): 2,  # interior InterVentricular-WM
            (5, 5, 2): 3,  # interior FornixWM
            (2, 6, 1): 4,  # lowercase alias: deliberately not selected
            (1, 4, 8): 5,  # missing hyphens: deliberately not selected
            (6, 0, 3): 6,  # extra hyphen: deliberately not selected
        },
    )
    aliases_only = _case(
        (5, 6, 7),
        {(0, 0, 0): 4, (2, 3, 5): 5, (4, 1, 6): 6},
    )
    # Positive diagonal entries encode array axes R, A, S; anisotropic voxel
    # sizes and nonzero origin catch accidental assumptions about spacing.
    affine = np.asarray(
        [[1.1, 0, 0, 17.0], [0, 1.7, 0, -23.0], [0, 0, 2.3, 8.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    save_nifti_ras(nib.Nifti1Image(seeded, affine), HERE / "seeded_labels.nii.gz", dtype=np.uint8)
    save_nifti_ras(
        nib.Nifti1Image(aliases_only, affine), HERE / "aliases_only_labels.nii.gz", dtype=np.uint8
    )


if __name__ == "__main__":
    main()
