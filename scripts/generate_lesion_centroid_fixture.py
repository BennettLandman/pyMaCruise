"""Create a compact synthetic RAS NIfTI for constrained centroid validation."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_constrained_centroids"


def main() -> None:
    # Frames: class labels; four memberships; two image channels; two fields.
    packed = np.zeros((5, 4, 3, 9), dtype=np.float32)
    ROOT.mkdir(parents=True, exist_ok=True)
    x, y, z = np.indices(packed.shape[:3])
    labels = ((x + 2 * y + z) % 3).astype(np.int16)
    packed[..., 0] = labels
    for k in range(4):
        steps = ((x * 3 + y * 5 + z * 7 + k * 2) % 9).astype(np.float32)
        packed[..., 1 + k] = np.float32(0.1) + steps * np.float32(0.05)
    for c in range(2):
        image_steps = ((x * 11 + y * 7 + z * 3 + c * 5) % 23).astype(np.float32)
        field_steps = ((x * 2 + y * 3 + z * 4 + c) % 5).astype(np.float32)
        packed[..., 5 + c] = np.float32(0.2) + image_steps * np.float32(0.1)
        packed[..., 7 + c] = np.float32(0.8) + field_steps * np.float32(0.1)
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
