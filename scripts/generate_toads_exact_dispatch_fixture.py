"""Generate a deterministic, asymmetric folded-shell RAS NIfTI test input."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_exact_dispatch"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nx, ny, nz = 15, 13, 11
    x, y, z = np.ogrid[:nx, :ny, :nz]
    dx = (x - 7.2) / 1.05
    dy = (y - 5.35) / 0.95
    dz = (z - 4.45) / 1.1
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    angle = np.arctan2(dy, dx)
    folded_inner = 3.2 + 0.65 * np.sin(5.0 * angle + 0.45 * dz)
    folded_outer = 5.55 + 0.72 * np.sin(5.0 * angle + 0.35 * dz + 0.4)
    ribbon = ((radius >= folded_inner) & (radius <= folded_outer)).astype(np.uint8)

    # Two channels carry deterministic intensity variation and real zeros for
    # the MASK mismatch branch. The shell has a left-right asymmetry and five
    # angular folds, so this is a geometry fixture rather than a uniform cube.
    image = np.zeros((nx, ny, nz, 2), dtype=np.float32)
    image[..., 0] = np.where(
        ribbon != 0,
        0.45 + 0.025 * x + 0.013 * y + 0.009 * z,
        np.where((x + 2 * y + z) % 7 == 0, 0.0, 0.08 + 0.01 * z),
    )
    image[..., 1] = np.where(
        ribbon != 0,
        0.72 - 0.018 * x + 0.011 * y + 0.006 * z,
        np.where((2 * x + y + 3 * z) % 9 == 0, 0.0, 0.14 + 0.008 * x),
    )
    affine = np.asarray(
        [[1.1, 0.0, 0.0, -7.7], [0.0, 1.2, 0.0, -6.0], [0.0, 0.0, 1.3, -5.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    header = nib.Nifti1Header()
    header.set_xyzt_units("mm")
    image_nii = nib.Nifti1Image(image, affine, header=header)
    image_nii.set_qform(affine, code=1)
    image_nii.set_sform(affine, code=1)
    nib.save(image_nii, OUTPUT / "images.nii.gz")
    ribbon_nii = nib.Nifti1Image(ribbon, affine)
    ribbon_nii.set_qform(affine, code=1)
    ribbon_nii.set_sform(affine, code=1)
    nib.save(ribbon_nii, OUTPUT / "folded_ribbon.nii.gz")
    print(f"wrote fixture to {OUTPUT}")


if __name__ == "__main__":
    main()
