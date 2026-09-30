"""Generate a synthetic RAS image for selected LesionTOADS variance formulas."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_variances"


def main() -> None:
    # Frames: labels, memberships, images, fields, lesion map, black-hole map.
    packed = np.zeros((5, 4, 3, 11), dtype=np.float32)
    ROOT.mkdir(parents=True, exist_ok=True)
    x, y, z = np.indices(packed.shape[:3])
    packed[..., 0] = ((x + 2 * y + z) % 3).astype(np.float32)
    for k in range(4):
        steps = ((x * 3 + y * 5 + z * 7 + k * 2) % 9).astype(np.float32)
        packed[..., 1 + k] = np.float32(0.1) + steps * np.float32(0.05)
    for c in range(2):
        image_steps = ((x * 11 + y * 7 + z * 3 + c * 5) % 23).astype(np.float32)
        field_steps = ((x * 2 + y * 3 + z * 4 + c) % 5).astype(np.float32)
        packed[..., 5 + c] = np.float32(0.2) + image_steps * np.float32(0.1)
        packed[..., 7 + c] = np.float32(0.8) + field_steps * np.float32(0.1)
    lesion_steps = ((x + 2 * y + z) % 9).astype(np.float32)
    black_steps = ((2 * x + y + 3 * z) % 9).astype(np.float32)
    lesion = np.float32(0.05) + lesion_steps * np.float32(0.06)
    black_hole = np.float32(0.04) + black_steps * np.float32(0.07)
    black_hole = np.where(((x + y + z) % 7) == 0, lesion, black_hole).astype(np.float32)
    zero = ((x + y + z) % 13) == 0
    packed[..., 9] = np.where(zero, 0.0, lesion)
    packed[..., 10] = np.where(zero, 0.0, black_hole)
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "variance_inputs.nii.gz"))


if __name__ == "__main__":
    main()
