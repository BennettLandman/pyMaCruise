"""Generate deterministic memberships for the TOADS smoothing slice."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_smoothing"
SHAPE = (7, 6, 5, 4)


def memberships() -> np.ndarray:
    x, y, z, class_index = np.indices(SHAPE, dtype=np.int32)
    return (((7 * x + 11 * y + 13 * z + 5 * class_index) % 47) / 46.0).astype(np.float32)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    image = nib.Nifti1Image(memberships(), affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, ROOT / "memberships.nii.gz")


if __name__ == "__main__":
    main()
