"""Generate small deterministic RAS NIfTI inputs for optimized TOADS likelihoods."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_optimized_likelihood"
SHAPE = (7, 6, 5, 3)


def inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, z, channel = np.indices(SHAPE, dtype=np.int32)
    images = (((11 * x + 7 * y + 5 * z + 3 * channel) % 101) / 100.0).astype(np.float32)
    field_mod = ((2 * x + 3 * y + 5 * z + channel) % 16).astype(np.float32)
    field = np.asarray(np.float32(0.85) + field_mod / np.float32(100.0), dtype=np.float32)
    bias_mod = ((5 * x + 3 * y + z + 2 * channel) % 17).astype(np.float32)
    bias = np.asarray(bias_mod / np.float32(1000.0), dtype=np.float32)
    return images, field, bias


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    for name, array in zip(("images", "field", "psf_bias"), inputs(), strict=True):
        image = nib.Nifti1Image(array, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
