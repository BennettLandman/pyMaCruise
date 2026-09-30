"""Generate a deterministic multi-channel RAS phantom for covariance tests."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_covariance"


def main() -> None:
    x, y, z = np.indices((5, 4, 3), dtype=np.int32)
    channels, classes = 3, 6
    packed = np.zeros((5, 4, 3, 1 + classes + channels + channels + 1), dtype=np.float32)
    labels = (x + 2 * y + z) % 5
    packed[..., 0] = labels.astype(np.float32)
    for k in range(classes):
        code = (x * 3 + y * 5 + z * 7 + k * 2) % 9
        membership = (code + 1).astype(np.float32) * np.float32(0.0625)
        if k < 5:
            membership = np.where((labels == k) & (((x + 2 * y + z) % 11) == 0), 0.0, membership)
        packed[..., 1 + k] = membership.astype(np.float32)
    for c in range(channels):
        image = ((x * 11 + y * 7 + z * 3 + c * 5) % 31).astype(np.float32) * np.float32(0.0625)
        field = np.float32(1.0) + ((x * 2 + y * 3 + z * 4 + c) % 5).astype(np.float32) * np.float32(0.0625)
        packed[..., 1 + classes + c] = image
        packed[..., 1 + classes + channels + c] = field
    lesion = (((x + 2 * y + z) % 7) + 1).astype(np.float32) * np.float32(0.0625)
    selected_tie = (labels == 2) & (((x + y + z) % 3) == 0)
    lesion = np.where(selected_tie, packed[..., 1 + 2], lesion).astype(np.float32)
    packed[..., -1] = lesion

    ROOT.mkdir(parents=True, exist_ok=True)
    affine = np.eye(4, dtype=np.float64)
    image = nib.Nifti1Image(packed, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, str(ROOT / "covariance_inputs.nii.gz"))


if __name__ == "__main__":
    main()
