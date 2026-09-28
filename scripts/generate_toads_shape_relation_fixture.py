"""Generate deterministic synthetic TOADS shape and relation NIfTIs."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_shape_relations"
SHAPE = (7, 6, 5, 4)


def inputs() -> tuple[np.ndarray, np.ndarray]:
    x, y, z, class_index = np.indices(SHAPE, dtype=np.int32)
    shape = (((x + 2 * y + 3 * z + 4 * class_index) % 11) / 10.0).astype(np.float32)
    labels = np.full(SHAPE, -1, dtype=np.int16)
    labels[..., 0] = ((x[..., 0] + y[..., 0] + z[..., 0]) % 4).astype(np.int16)
    direct2 = ((x[..., 0] + y[..., 0] + z[..., 0]) % 3) == 0
    labels[..., 1][direct2] = (labels[..., 0][direct2] + 1) % 4
    labels[..., 2] = ((2 * x[..., 0] + y[..., 0]) % 4).astype(np.int16)
    labels[..., 3] = ((x[..., 0] + 2 * y[..., 0] + z[..., 0] + 1) % 4).astype(np.int16)
    return shape, labels


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    for name, array in zip(("shape_priors", "closest_labels"), inputs(), strict=True):
        image = nib.Nifti1Image(array, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        nib.save(image, ROOT / f"{name}.nii.gz")


if __name__ == "__main__":
    main()
