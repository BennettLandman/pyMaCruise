"""Generate the deterministic folded-cortex-like TOADS shape prior fixture."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests" / "data" / "synthetic" / "toads_shape_bounds_folded"
SHAPE = (39, 35, 31)


def folded_prior() -> np.ndarray:
    x, y, z = np.indices(SHAPE, dtype=np.float64)
    dx = x - (SHAPE[0] - 1) / 2.0
    dy = y - (SHAPE[1] - 1) / 2.0
    dz = z - (SHAPE[2] - 1) / 2.0
    radius = np.sqrt(dx * dx + dy * dy + 0.35 * dz * dz)
    fold = (
        2.4 * np.sin(dx * 0.48 + dz * 0.19)
        + 1.7 * np.sin(dy * 0.53 - dz * 0.31)
        + 1.1 * np.sin((dx + dy) * 0.29 + dz * 0.37)
    )
    distance = np.abs(radius - (12.0 + fold))
    return np.asarray(np.maximum(0.0, 1.0 - distance / 2.6), dtype=np.float32)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(folded_prior(), np.diag([1.0, 1.0, 1.0, 1.0]))
    image.set_qform(image.affine, code=1)
    image.set_sform(image.affine, code=1)
    nib.save(image, ROOT / "folded_shape_prior.nii.gz")


if __name__ == "__main__":
    main()
