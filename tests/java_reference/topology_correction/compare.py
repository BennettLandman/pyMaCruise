"""Report exact-source Java vs Python topology-correction differences."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.topology import topology_correct

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "synthetic" / "topology_correction"
FIXTURES = (
    "spatial_hash",
    "plateau",
    "boundary",
    "narrow_branch",
    "disconnected_seeds",
)


def main() -> None:
    total_voxels = 0
    total_different = 0
    for fixture in FIXTURES:
        fixture_dir = DATA_DIR if fixture == "spatial_hash" else DATA_DIR / fixture
        source = nib.load(str(fixture_dir / "input.nii.gz"))
        paint = nib.load(str(fixture_dir / "paint_mask.nii.gz"))
        java = np.asarray(
            nib.load(str(fixture_dir / "java_reference.nii.gz")).dataobj,
            dtype=np.float32,
        )
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
        total_voxels += absolute.size
        total_different += int(np.count_nonzero(absolute))
        print(
            f"{fixture}: voxels={absolute.size} "
            f"different={np.count_nonzero(absolute)} "
            f"max_absolute_difference={absolute.max():.9g} "
            f"java_sum={java.sum(dtype=np.float64):.9g} "
            f"python_sum={python.sum(dtype=np.float64):.9g}"
        )
    print(f"total_voxels={total_voxels}")
    print(f"total_different_voxels={total_different}")


if __name__ == "__main__":
    main()
