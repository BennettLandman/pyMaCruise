"""Create the deterministic, asymmetric default-TOADS RAS NIfTI phantom."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
AFFINE = np.diag([0.9, 1.2, 1.45, 1.0])


def main() -> None:
    x, y, z = np.indices((7, 6, 5), dtype=np.float32)
    images = np.empty((7, 6, 5, 2), dtype=np.float32)
    images[..., 0] = 0.17 + 0.013 * x + 0.007 * y + 0.009 * z + 0.001 * x * z
    images[..., 1] = 0.29 + 0.004 * x + 0.011 * y + 0.006 * z + 0.001 * y * z

    shapes = np.empty((7, 6, 5, 11), dtype=np.float32)
    for klass in range(11):
        ridge = np.abs((x + 2 * y + 3 * z + klass) % 9 - 4)
        band = np.float32(0.97 - 0.13 * ridge)
        class_gain = np.float32(0.08 + 0.027 * klass)
        mask = ((x + 2 * y + z + klass) % 5) != 0
        shapes[..., klass] = np.where(mask, np.maximum(np.float32(0.0), band * class_gain), 0.0)

    for name, data in (("images", images), ("transformed_shapes", shapes)):
        image = nib.Nifti1Image(data, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
