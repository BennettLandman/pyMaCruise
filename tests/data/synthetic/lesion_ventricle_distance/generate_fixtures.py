"""Create the asymmetric true-RAS NIfTI fixture for selected ventricle names."""

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import save_nifti_ras


HERE = Path(__file__).resolve().parent


def main() -> None:
    seeded = np.zeros((7, 8, 9), dtype=np.uint8)
    for point, label in {
        (0, 2, 7): 1,  # border Ventricles
        (4, 1, 3): 2,  # interior LeftVentricle
        (5, 5, 2): 3,  # interior RightVentricle
        (6, 7, 0): 4,  # border ThirdVentricle
        (2, 6, 1): 5,  # lowercase negative control
        (1, 4, 8): 6,  # punctuation negative control
        (3, 0, 5): 7,  # mixed-case negative control
        (6, 0, 3): 8,  # unrelated negative control
    }.items():
        seeded[point] = label

    nonmatching = np.zeros((5, 6, 7), dtype=np.uint8)
    for point, label in {
        (0, 0, 0): 5,
        (2, 3, 5): 6,
        (4, 1, 6): 7,
        (3, 2, 1): 8,
    }.items():
        nonmatching[point] = label

    # Positive diagonal entries encode true RAS voxel axes; nonzero origin and
    # anisotropic spacing make the fixture's spatial header unambiguous.
    affine = np.asarray(
        [[1.1, 0, 0, 17.0], [0, 1.7, 0, -23.0], [0, 0, 2.3, 8.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    save_nifti_ras(nib.Nifti1Image(seeded, affine), HERE / "seeded_labels.nii.gz", dtype=np.uint8)
    save_nifti_ras(
        nib.Nifti1Image(nonmatching, affine),
        HERE / "nonmatching_labels.nii.gz",
        dtype=np.uint8,
    )


if __name__ == "__main__":
    main()
