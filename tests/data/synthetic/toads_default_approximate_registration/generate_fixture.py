"""Generate asymmetric folded RAS volumes for the default TOADS pyramid."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
AFFINE = np.diag([.9, 1.2, 1.45, 1.0])


def folded_shapes(x: np.ndarray, y: np.ndarray, z: np.ndarray, phase: int) -> np.ndarray:
    maps = np.empty(x.shape + (11,), dtype=np.float32)
    for klass in range(11):
        fold = np.abs((x + 2 * y + 3 * z + 2 * klass + phase) % 11 - 5)
        amplitude = np.maximum(np.float32(0.0), np.float32(.94) - np.float32(.11) * fold)
        maps[..., klass] = amplitude * np.float32(.09 + .022 * klass)
    return maps


def save(name: str, values: np.ndarray) -> None:
    image = nib.Nifti1Image(values, AFFINE)
    image.set_qform(AFFINE, code=1)
    image.set_sform(AFFINE, code=1)
    nib.save(image, ROOT / f"{name}.nii.gz")


def main() -> None:
    x, y, z = np.indices((9, 8, 7), dtype=np.float32)
    images = np.empty((9, 8, 7, 2), dtype=np.float32)
    images[..., 0] = .14 + .009 * x + .006 * y + .008 * z + .0007 * x * z
    images[..., 1] = .31 + .003 * x + .008 * y + .005 * z + .0006 * y * z
    save("images", images)
    save("transformed_shapes", folded_shapes(x, y, z, phase=1))

    tx, ty, tz = np.indices((11, 10, 9), dtype=np.float32)
    save("atlas_shapes", folded_shapes(tx, ty, tz, phase=3))


if __name__ == "__main__":
    main()
