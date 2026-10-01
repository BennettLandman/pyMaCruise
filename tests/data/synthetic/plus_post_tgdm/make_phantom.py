from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
SOURCE = HERE.parent / "plus_surfaces"


def main() -> None:
    """Materialize fixed 13^3 cropped post-TGDM inputs from the RAS level sets."""
    names = ("initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer")
    for name in names:
        if name in {"initial_level_set", "outer"}:
            source = nib.load(SOURCE / "outer.nii.gz")
            data = np.asarray(source.dataobj)
        elif name in {"initial_inner", "inner"}:
            source = nib.load(SOURCE / "inner.nii.gz")
            data = np.asarray(source.dataobj)
        elif name == "central":
            source = nib.load(SOURCE / "central.nii.gz")
            data = np.asarray(source.dataobj)
        else:
            source = nib.load(SOURCE / "inner.nii.gz")
            data = np.zeros(source.shape, dtype=np.uint8)
        cropped = data[2:15, 2:15, 2:15]
        affine = source.affine.copy()
        translation = np.eye(4)
        translation[:3, 3] = (2, 2, 2)
        affine = affine @ translation
        image = nib.Nifti1Image(cropped, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, HERE / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
