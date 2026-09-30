"""Generate a compact synthetic RAS input for LesionTOADS channel weights."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_channel_weights"


def main() -> None:
    # Frames: six class memberships, two image channels, two fields, and three
    # auxiliary maps (lesion, black-hole, outlier membership).
    packed = np.zeros((5, 4, 3, 13), dtype=np.float32)
    ROOT.mkdir(parents=True, exist_ok=True)
    x, y, z = np.indices(packed.shape[:3])
    for k in range(6):
        steps = ((x * 3 + y * 5 + z * 7 + k * 2) % 9).astype(np.float32)
        packed[..., k] = np.float32(0.1) + steps * np.float32(0.05)
    for c in range(2):
        image_steps = ((x * 11 + y * 7 + z * 3 + c * 5) % 23).astype(np.float32)
        field_steps = ((x * 2 + y * 3 + z * 4 + c) % 5).astype(np.float32)
        packed[..., 6 + c] = np.float32(0.2) + image_steps * np.float32(0.1)
        packed[..., 8 + c] = np.float32(0.8) + field_steps * np.float32(0.1)
    packed[..., 10] = np.float32(0.1) + ((x + 2 * y + z) % 7).astype(np.float32) * np.float32(0.05)
    packed[..., 11] = np.float32(0.08) + ((2 * x + y + 3 * z) % 8).astype(np.float32) * np.float32(0.04)
    packed[..., 12] = np.float32(0.05) + ((x + 3 * y + 2 * z) % 7).astype(np.float32) * np.float32(0.06)
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "channel_weight_inputs.nii.gz"))


if __name__ == "__main__":
    main()
