"""Generate deterministic RAS NIfTIs for the TOADS PSF ridge-bias slice."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_psf_bias"
SHAPE = (9, 8, 7, 2)


def inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, z, channel = np.indices(SHAPE, dtype=np.int32)
    images = (((13 * x + 11 * y + 7 * z + 5 * channel) % 113) / 112.0).astype(np.float32)
    field_remainder = ((3 * x + 5 * y + 2 * z + channel) % 21).astype(np.float32)
    field = np.asarray(np.float32(0.8) + field_remainder / np.float32(100.0), dtype=np.float32)
    labels = np.broadcast_to((x >= 4).astype(np.int16), SHAPE).copy()[..., 0]
    labels[2, 3, 3] = 2
    labels[6, 1:4, 2:5] = 0
    return images, field, labels


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    images, field, labels = inputs()
    for name, array in (("images", images), ("field", field), ("segmentation", labels)):
        image = nib.Nifti1Image(array, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
