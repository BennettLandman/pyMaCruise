"""Regenerate the committed 17-cubed synthetic Dura inputs as RAS NIfTI."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (17, 17, 17)
AFFINE_RAS = np.eye(4, dtype=np.float64)


def main() -> None:
    x, y, z = np.indices(SHAPE, dtype=np.float64)
    radius = np.sqrt((x - 8.0) ** 2 + (y - 8.0) ** 2 + (z - 8.0) ** 2)
    white_matter = (255.0 * np.clip((4.4 - radius) / 1.7, 0.0, 1.0)).astype(
        np.float32
    )
    gray_matter = (
        255.0 * np.clip(1.0 - np.abs(radius - 4.25) / 1.8, 0.0, 1.0)
    ).astype(np.float32)
    t1 = (80.0 + 2.0 * x + 3.0 * y + 5.0 * z).astype(np.float32)
    brain_mask = (radius <= 2.9).astype(np.uint8)

    for case_name in ("attached_dura", "no_dura"):
        case_dir = ROOT / case_name
        case_dir.mkdir(parents=True, exist_ok=True)
        case_gray_matter = gray_matter.copy()
        if case_name == "attached_dura":
            case_gray_matter[(y == 8) & (z == 8) & (x >= 12) & (x <= 15)] = 190.0
            case_gray_matter[
                (x >= 14) & (x <= 15) & (y >= 7) & (y <= 9) & (z >= 7) & (z <= 9)
            ] = 230.0
        for name, data in (
            ("t1", t1),
            ("wm", white_matter),
            ("gm", case_gray_matter),
            ("brain_mask", brain_mask),
        ):
            image = nib.Nifti1Image(data, AFFINE_RAS)
            image.set_qform(AFFINE_RAS, code=1)
            image.set_sform(AFFINE_RAS, code=1)
            nib.save(image, str(case_dir / f"{name}.nii.gz"))


if __name__ == "__main__":
    main()
