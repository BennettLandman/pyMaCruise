"""Generate a small RAS transformed-shape volume for TOADS prior tests."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_inverse_shape_prior"
SHAPE = (5, 6, 4, 3)


def main() -> None:
    x, y, z, class_index = np.indices(SHAPE, dtype=np.int32)
    shape_values = np.asarray(((3 * x + 5 * y + 7 * z + 11 * class_index) % 23) / 22.0, dtype=np.float32)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    image = nib.Nifti1Image(shape_values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    ROOT.mkdir(parents=True, exist_ok=True)
    nib.save(image, ROOT / "shape_values.nii.gz")


if __name__ == "__main__":
    main()
