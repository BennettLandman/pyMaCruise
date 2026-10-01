from __future__ import annotations

import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
SHAPE = (24, 24, 24)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def main() -> None:
    grid = np.indices(SHAPE, dtype=np.float32)
    radius = np.sqrt(sum((axis - 11.5) ** 2 for axis in grid))
    white = np.where(radius <= 4.0, 230.0, 0.0).astype(np.float32)
    gray = np.where((radius > 5.0) & (radius <= 9.0), 220.0, 0.0).astype(np.float32)
    brain = (radius <= 10.0).astype(np.uint8)

    label200_control = np.zeros(SHAPE, dtype=np.int16)
    label200_control[(radius <= 4.0)] = 201
    label200_control[(radius > 4.0) & (radius <= 10.0)] = 200
    label202_overlook_patch = label200_control.copy()
    # The patch changes to selected boundary-overlook label 202; this is not
    # the LesionTOADS atlas class named type=lesion.
    label202_overlook_patch[8:12, 14:18, 10:14] = 202
    toads_wm_boundary = np.where(radius <= 4.0, 0.5, 0.0).astype(np.float32)

    files = {
        "segmentation_label202_overlook_patch.nii.gz": label202_overlook_patch,
        "segmentation_label200_control.nii.gz": label200_control,
        "toads_wm_boundary.nii.gz": toads_wm_boundary,
        "multi_atlas_gm.nii.gz": gray.astype(np.uint8),
        "multi_atlas_wm.nii.gz": white.astype(np.uint8),
        "multi_atlas_brain_mask.nii.gz": brain,
        "toads_gm.nii.gz": np.where(gray > 0, 0.1, 0.0).astype(np.float32),
        "toads_wm.nii.gz": np.where(white > 0, 0.6, 0.0).astype(np.float32),
        "toads_mask.nii.gz": np.where(brain, 7, 0).astype(np.int16),
    }
    for name, data in files.items():
        image = nib.Nifti1Image(data, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, HERE / name)

    names = sorted([*files, "README.md", "make_phantom.py"])
    lines = []
    for name in names:
        lines.append(
            f"{hashlib.sha256((HERE / name).read_bytes()).hexdigest()}  {name}"
        )
    (HERE / "SHA256SUMS").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
