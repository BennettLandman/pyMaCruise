"""Convert selected Java two-call Dura topology outputs to NIfTI fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = CASE_DIR.parents[1] / "data" / "synthetic" / "dura_topology_two_pass"
SHAPE = (9, 9, 9)
SECTIONS = {
    "WM_INPUT_BEGIN": "wm_input.nii.gz",
    "PAINT_BEGIN": "paint_mask.nii.gz",
    "RAW_GM_BEGIN": "raw_gm.nii.gz",
    "WM_TC_BEGIN": "java_wm_tc.nii.gz",
    "UNION_BEGIN": "java_union.nii.gz",
    "GM_TC_BEGIN": "java_gm_tc.nii.gz",
}


def main() -> None:
    lines = Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
    headings = list(SECTIONS)
    affine = np.eye(4, dtype=np.float64)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for index, heading in enumerate(headings):
        start = lines.index(heading) + 1
        next_start = next(
            (lines.index(candidate) for candidate in headings[index + 1 :]), len(lines)
        )
        values = np.asarray([float(value) for value in lines[start:next_start]], dtype=np.float32)
        if values.size != int(np.prod(SHAPE)):
            raise ValueError(f"{heading} has {values.size} voxels; expected {np.prod(SHAPE)}")
        # Java ModelImage is x-fast; restore that order for NIfTI x/y/z axes.
        data = values.reshape(SHAPE, order="F")
        nib.save(nib.Nifti1Image(data, affine), OUTPUT_DIR / SECTIONS[heading])


if __name__ == "__main__":
    main()
