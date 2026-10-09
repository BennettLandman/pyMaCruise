#!/usr/bin/env python3
"""Create the committed folded cortical label phantom for boundary comparison."""
from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    nx, ny, nz = 66, 76, 66
    x, y, z = np.meshgrid(
        np.linspace(-1, 1, nx, dtype=np.float32),
        np.linspace(-1, 1, ny, dtype=np.float32),
        np.linspace(-1, 1, nz, dtype=np.float32),
        indexing="ij",
    )
    radius = np.sqrt((x / 0.78) ** 2 + (y / 0.86) ** 2 + (z / 0.82) ** 2)
    labels = np.zeros((nx, ny, nz), dtype=np.uint16)
    labels[radius < 0.69] = 110
    cortex = (radius >= 0.69) & (radius <= 0.92)
    angle = np.arctan2(y, x)
    folded = np.sin(7 * angle + 2.5 * z + 1.8 * np.sin(3 * angle - z))
    labels[cortex & (folded >= 0)] = 100
    labels[cortex & (folded < 0)] = 101
    low_region = ((x + 0.35) ** 2 + (y - 0.04) ** 2 + (z + 0.08) ** 2 < 0.075) & cortex
    labels[low_region] = 4
    overlook = ((x - 0.58) ** 2 + (y + 0.12) ** 2 + z**2 < 0.025) & cortex
    labels[overlook] = 120
    affine = np.array([
        [0.3047619163990021, 0, 0, -0.3047619163990021],
        [0, 0.24303798377513885, 0, 0.9721519351005554],
        [0, 0, 0.4571428894996643, -0.4571428894996643],
        [0, 0, 0, 1],
    ])
    out = nib.Nifti1Image(labels, affine)
    out.set_data_dtype(np.uint16)
    out.header.set_slope_inter(None, None)
    destination = root / "tests/data/synthetic/plus_boundary_folded_cortex/segmentation.nii.gz"
    destination.parent.mkdir(parents=True, exist_ok=True)
    nib.save(out, destination)
    print(f"wrote {destination} shape={labels.shape} labels={dict(zip(*np.unique(labels, return_counts=True), strict=True))}")


if __name__ == "__main__":
    main()
