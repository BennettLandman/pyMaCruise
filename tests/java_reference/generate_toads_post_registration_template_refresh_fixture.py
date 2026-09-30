"""Create the synthetic RAS NIfTI pair for the P18czr method harness."""

from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).parents[1] / "data" / "synthetic" / "toads_post_registration_template_refresh"


def image_value(c: int, x: int, y: int, z: int) -> np.float32:
    dx, dy, dz = x - 6, y - 6, z - 5
    if dx * dx + dy * dy + dz * dz > 29 or (x + 2 * y + 3 * z + c) % 7 == 0:
        return np.float32(0.0)
    return np.float32(((x * 17 + y * 11 + z * 5 + c * 13) % 29 + 1) / 32.0)


def atlas_label(x: int, y: int, z: int) -> np.uint8:
    shell = (
        3 <= x <= 12 and 2 <= y <= 13 and 2 <= z <= 12
        and (abs(y - (5 + x % 4)) <= 1 or abs(z - (4 + x % 5)) <= 1)
    )
    if shell:
        return np.uint8(2)
    if (
        6 <= x <= 14 and 7 <= y <= 14 and 5 <= z <= 13
        and ((x + y + z) % 4 != 0 or abs(y - (8 + x % 3)) <= 1)
    ):
        return np.uint8(3)
    if x == 2 and 8 <= y <= 11 and z == 10:
        return np.uint8(4)
    return np.uint8(1)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    images = np.empty((14, 13, 11, 2), dtype=np.float32)
    for x in range(14):
        for y in range(13):
            for z in range(11):
                for channel in range(2):
                    images[x, y, z, channel] = image_value(channel, x, y, z)
    topology = np.empty((17, 16, 15), dtype=np.uint8)
    for x in range(17):
        for y in range(16):
            for z in range(15):
                topology[x, y, z] = atlas_label(x, y, z)
    nib.save(
        nib.Nifti1Image(images, np.diag([0.8, 1.1, 1.4, 1.0])),
        ROOT / "native_t1_channels_ras.nii.gz",
    )
    nib.save(
        nib.Nifti1Image(topology, np.diag([1.0, 1.2, 1.5, 1.0])),
        ROOT / "atlas_template_ras.nii.gz",
    )


if __name__ == "__main__":
    main()
