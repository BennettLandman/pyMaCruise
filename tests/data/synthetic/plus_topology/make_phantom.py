from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
SOURCE = HERE.parent / "topology_correction"


def main() -> None:
    volume = nib.load(SOURCE / "input.nii.gz")
    paint = nib.load(SOURCE / "paint_mask.nii.gz")
    nib.save(volume, HERE / "topology_volume.nii.gz")
    nib.save(paint, HERE / "paint_mask.nii.gz")
    ma_wm = np.zeros(volume.shape, dtype=np.uint8)
    ma_wm[2:7, 2:7, 2:7] = 255
    wm = nib.Nifti1Image(ma_wm, volume.affine)
    wm.set_qform(volume.affine, code=1)
    wm.set_sform(volume.affine, code=1)
    nib.save(wm, HERE / "ma_wm.nii.gz")


if __name__ == "__main__":
    main()
