"""Generate the categorical RAS fixture for optimized TOADS class variances."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_optimized_class_variance"
SHAPE = (7, 6, 5)


def main() -> None:
    x, y, z = np.indices(SHAPE, dtype=np.int16)
    labels = np.asarray((x + 2 * y + z) % 3, dtype=np.int16)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    image = nib.Nifti1Image(labels, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    ROOT.mkdir(parents=True, exist_ok=True)
    nib.save(image, ROOT / "segmentation.nii.gz")


if __name__ == "__main__":
    main()
