"""Convert the selected Java topology source's stdout into local NIfTI fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = CASE_DIR.parents[1] / "data" / "synthetic" / "topology_correction"
SHAPE = (9, 9, 9)


def main() -> None:
    lines = Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
    values = np.asarray([float(line) for line in lines[lines.index("OUTPUT_BEGIN") + 1 :]], dtype=np.float32)
    if values.size != int(np.prod(SHAPE)):
        raise ValueError(f"Java output has {values.size} values; expected {np.prod(SHAPE)}")

    x, y, z = np.indices(SHAPE, dtype=np.int64)
    hash_values = ((x * 73856093) ^ (y * 19349663) ^ (z * 83492791)) & 0xFFFFFFFF
    intensity = ((hash_values & 0x7FFFFFFF) % 100000).astype(np.float32) / 100000.0
    intensity[0, 0, 0] = 0.0
    paint = np.zeros(SHAPE, dtype=np.float32)
    paint[4, 4, 4] = 1.0
    paint[4, 4, 5] = 0.5

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    affine = np.eye(4, dtype=np.float64)
    for name, data in (
        ("input.nii.gz", intensity),
        ("paint_mask.nii.gz", paint),
        ("java_reference.nii.gz", values.reshape(SHAPE)),
    ):
        nib.save(nib.Nifti1Image(data, affine), OUTPUT_DIR / name)


if __name__ == "__main__":
    main()
