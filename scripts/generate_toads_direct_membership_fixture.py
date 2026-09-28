"""Generate RAS NIfTI energy arrays for direct TOADS membership updates."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_direct_membership"
SHAPE = (7, 6, 5, 4)


def inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, z, class_index = np.indices(SHAPE, dtype=np.int32)
    energy = (((13 * x + 7 * y + 3 * z + 11 * class_index) % 29) + 1).astype(np.float32) / np.float32(13.0)
    objective = (((5 * x + 3 * y + 7 * z + 2 * class_index) % 19) + 1).astype(np.float32) / np.float32(17.0)
    initial = (((5 * x + 7 * y + 11 * z + 3 * class_index) % 31) / 30.0).astype(np.float32)
    threshold = np.float32(1e-10)
    energy[2, 2, 2, 0] = np.float32(0.0)
    energy[3, 2, 2, 1] = threshold
    energy[4, 3, 2, 2] = np.nextafter(threshold, np.float32(0.0))
    energy[2, 3, 2, 3] = np.nextafter(threshold, np.float32(1.0))
    return energy, objective, initial


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    for name, array in zip(("data_energy", "objective_energy", "initial_memberships"), inputs(), strict=True):
        image = nib.Nifti1Image(array, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
