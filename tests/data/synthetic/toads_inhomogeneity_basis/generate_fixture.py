"""Generate the small RAS NIfTI grid used by the Chebyshev-basis harness."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
AFFINE = np.diag([.8, 1.1, 1.7, 1.0])


def main() -> None:
    x, y, z = np.indices((5, 4, 3), dtype=np.float32)
    data = np.float32(.2) + np.float32(.03) * x + np.float32(.04) * y + np.float32(.05) * z
    image = nib.Nifti1Image(data, AFFINE)
    image.set_qform(AFFINE, code=1)
    image.set_sform(AFFINE, code=1)
    nib.save(image, ROOT / "grid.nii.gz")


if __name__ == "__main__":
    main()
