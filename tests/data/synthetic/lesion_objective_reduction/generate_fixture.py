"""Generate the synthetic NIfTI inputs for selected-source objective reduction."""

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
SHAPE = (5, 6, 7)


def make_fixture() -> np.ndarray:
    volume = np.zeros((*SHAPE, 11), dtype=np.float32)
    volume[..., :3] = np.asarray([0.5, 0.25, 0.25], dtype=np.float32)
    volume[..., 3:6] = np.asarray([64.0, 16.0, 16.0], dtype=np.float32)
    volume[..., 6] = np.float32(1.0)
    volume[..., 7] = np.float32(0.25)
    volume[..., 8] = np.float32(0.125)
    volume[..., 9] = np.float32(2.0)
    volume[..., 10] = np.float32(4.0)

    # Zero class terms still leave the auxiliary objective terms active.
    volume[1, 1, 1, 3:6] = np.float32(0.0)
    # A zero denominator suppresses classes, but source adds auxiliary terms.
    volume[1, 1, 2, 6] = np.float32(0.0)
    volume[1, 1, 2, 7] = np.float32(0.5)
    volume[1, 1, 2, 8] = np.float32(0.25)
    # A negative or NaN denominator also fails Java's strict `den > 0` check.
    volume[1, 1, 3, 6] = np.float32(-1.0)
    volume[1, 1, 4, 6] = np.float32(np.nan)
    # IEEE infinity followed by Inf*0 (NaN), then an explicit NaN class energy.
    volume[3, 4, 3, 5] = np.float32(np.inf)
    volume[3, 4, 4, 5] = np.float32(np.inf)
    volume[3, 4, 4, 2] = np.float32(0.0)
    volume[3, 4, 5, 4] = np.float32(np.nan)
    return volume


def main() -> None:
    affine = np.diag([1.1, 1.3, 1.7, 1.0]).astype(np.float64)
    image = nib.Nifti1Image(make_fixture(), affine)
    image.header.set_xyzt_units("mm")
    nib.save(image, HERE / "objective_inputs_ras.nii.gz")


if __name__ == "__main__":
    main()
