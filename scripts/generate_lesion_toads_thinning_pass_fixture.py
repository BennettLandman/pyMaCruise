"""Generate the selected one-pass LesionTOADS synthetic RAS input volume."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_thinning_pass"
N = 11


def main() -> None:
    # Frames: segmentation, initial skeleton, class memberships 0..2, lesion,
    # black hole, transformed shapes 0..2, initial availability 0..2, ordering.
    packed = np.zeros((N, N, N, 14), dtype=np.float32)
    x, y, z = np.indices((N, N, N))
    label = np.where(x <= 4, 0, np.where(x <= 7, 1, 2)).astype(np.int16)
    packed[..., 0] = label
    border = ((x == 0) | (x == N - 1) | (y == 0) | (y == N - 1) | (z == 0) | (z == N - 1))
    packed[..., 1] = np.where(border, 7, -1)
    for class_index in range(3):
        values = np.where(label == class_index, 0.8, 0.05).astype(np.float32)
        if class_index == 1:
            values = np.where(label == 1, 0.8, 0.1).astype(np.float32)
        packed[..., 2 + class_index] = values
        packed[..., 7 + class_index] = 0.0
        packed[..., 10 + class_index] = ((x + y + z + class_index) % 2 == 0).astype(np.float32)
    packed[..., 5] = 0.1
    packed[..., 6] = 0.05
    packed[..., 13] = 0.375
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "thinning_pass_inputs.nii.gz"))


if __name__ == "__main__":
    main()
