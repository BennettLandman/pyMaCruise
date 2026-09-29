"""Report exact-source Java vs Python topology-correction differences."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.topology import topology_correct

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "synthetic" / "topology_correction"


def main() -> None:
    load = lambda name: nib.load(str(DATA_DIR / name))
    source = load("input.nii.gz")
    paint = load("paint_mask.nii.gz")
    java = np.asarray(load("java_reference.nii.gz").dataobj, dtype=np.float32)
    python = np.asarray(
        topology_correct(
            source,
            paint,
            paint_threshold=0.5,
            lowest_level=0.0,
            highest_level=1.0,
            regularize=0.05,
        ).corrected.dataobj,
        dtype=np.float32,
    )
    absolute = np.abs(java - python)
    print(f"voxels={absolute.size}")
    print(f"different_voxels={np.count_nonzero(absolute)}")
    print(f"max_absolute_difference={absolute.max():.9g}")
    print(f"java_sum={java.sum(dtype=np.float64):.9g}")
    print(f"python_sum={python.sum(dtype=np.float64):.9g}")


if __name__ == "__main__":
    main()
