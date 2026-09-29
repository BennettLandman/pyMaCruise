from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent


def main() -> None:
    shape = (7, 7, 7)
    affine = np.diag([1.2, 1.3, 1.4, 1.0])
    atlas_gm = np.zeros(shape, dtype=np.uint8)
    atlas_wm = np.zeros(shape, dtype=np.uint8)
    brain_mask = np.ones(shape, dtype=np.uint8)
    toads_gm = np.full(shape, 0.5, dtype=np.float32)
    toads_wm = np.full(shape, 0.6, dtype=np.float32)
    toads_mask = np.full(shape, 7, dtype=np.int16)

    atlas_wm[1:6, 1:6, 1:6] = 255
    atlas_gm[1, 1, 1] = 255
    atlas_gm[1, 2, 2] = 255
    toads_gm[3, 3, 3] = 0.1  # confident atlas WM has priority
    toads_gm[1, 1, 1] = 0.1  # confident GM fill, outside eroded WM
    toads_gm[1, 2, 2] = 0.255  # strict threshold: fallback at equality
    brain_mask[0, 0, 0] = 0  # outside-brain zeroing

    arrays = {
        "multi_atlas_gm": atlas_gm,
        "multi_atlas_wm": atlas_wm,
        "multi_atlas_brain_mask": brain_mask,
        "toads_gm": toads_gm,
        "toads_wm": toads_wm,
        "toads_mask": toads_mask,
    }
    for name, data in arrays.items():
        image = nib.Nifti1Image(data, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, HERE / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
