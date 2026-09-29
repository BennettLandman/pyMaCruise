from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent


def main() -> None:
    shape = (17, 17, 17)
    affine = np.diag([1.0, 1.5, 2.0, 1.0])
    i, j, k = np.indices(shape, dtype=np.float32)
    distance_mm = np.sqrt((i - 8) ** 2 + ((j - 8) * 1.5) ** 2 + ((k - 8) * 2.0) ** 2)
    for name, radius_mm in (("inner", 3.0), ("central", 4.0), ("outer", 5.0)):
        image = nib.Nifti1Image((distance_mm - radius_mm).astype(np.float32), affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, HERE / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
