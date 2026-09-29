from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent


def main() -> None:
    shape = (13, 13, 13)
    affine = np.diag([1.0, 1.25, 1.5, 1.0])
    segmentation = np.zeros(shape, dtype=np.int16)
    segmentation[2:11, 2:11, 2:11] = 200
    segmentation[5:9, 5:9, 5:9] = 201
    segmentation[3:5, 4:9, 4:9] = 120  # selected overlooked label
    toads_wm = np.zeros(shape, dtype=np.float32)
    toads_wm[7:11, 2:7, 2:11] = 0.5  # cleanup removes boundary in this region
    for name, data in (("segmentation", segmentation), ("toads_wm_mask", toads_wm)):
        image = nib.Nifti1Image(data, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, HERE / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
