"""Convert selected Java source stdout into RAS NIfTI reference fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
DATA_DIR = CASE_DIR.parents[1] / "data" / "synthetic" / "topology_correction"
SHAPES = {
    "spatial_hash": (9, 9, 9),
    "plateau": (9, 9, 9),
    "boundary": (9, 9, 9),
    "narrow_branch": (11, 9, 9),
    "disconnected_seeds": (9, 9, 9),
}
SECTIONS = ("INPUT_BEGIN", "PAINT_BEGIN", "OUTPUT_BEGIN")


def main() -> None:
    fixture_name = sys.argv[1]
    output_path = Path(sys.argv[2])
    shape = SHAPES[fixture_name]
    lines = output_path.read_text(encoding="utf-8").splitlines()
    sections: dict[str, np.ndarray] = {}
    for index, heading in enumerate(SECTIONS):
        start = lines.index(heading) + 1
        next_heading = next(
            (lines.index(candidate) for candidate in SECTIONS[index + 1 :]), len(lines)
        )
        values = np.asarray([float(line) for line in lines[start:next_heading]], dtype=np.float32)
        if values.size != int(np.prod(shape)):
            raise ValueError(
                f"{fixture_name} {heading}: got {values.size} values, expected {np.prod(shape)}"
            )
        # Java ModelImage's x coordinate is the fastest varying index. Restore
        # that x/y/z array order rather than NumPy's default C-order reshape.
        sections[heading] = values.reshape(shape, order="F")

    destination = DATA_DIR if fixture_name == "spatial_hash" else DATA_DIR / fixture_name
    destination.mkdir(parents=True, exist_ok=True)
    affine = np.eye(4, dtype=np.float64)
    for heading, filename in zip(SECTIONS, ("input.nii.gz", "paint_mask.nii.gz", "java_reference.nii.gz"), strict=True):
        nib.save(nib.Nifti1Image(sections[heading], affine), destination / filename)


if __name__ == "__main__":
    main()
