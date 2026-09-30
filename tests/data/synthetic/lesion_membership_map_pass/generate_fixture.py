"""Generate the committed deterministic folded RAS pass input fixture."""

from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
SHAPE = (7, 8, 9)


def make_fixture() -> np.ndarray:
    x, y, z = np.indices(SHAPE, dtype=np.float32)
    channels: list[np.ndarray] = []
    memberships = []
    class_data = []
    shape_values = []
    for k in range(5):
        memberships.append(
            np.asarray(0.25 + x / 64 + y / 128 + z / 256 + k / 32, dtype=np.float32)
        )
        class_data.append(
            np.asarray(0.75 + x / 16 + y / 32 + z / 64 + k / 16, dtype=np.float32)
        )
        folded_shape = np.asarray(0.25 + x / 32 + y / 64 + z / 128 + k / 32, dtype=np.float32)
        fold_limit = 2 + ((2 * y.astype(np.int32) + z.astype(np.int32) + (y * z).astype(np.int32)) % 4)
        folded_shape = np.where(x <= fold_limit, folded_shape, np.float32(0.0)).astype(np.float32)
        shape_values.append(folded_shape)

    # Fixed sparse support holes also cover the source's shape<=0 clearing branch.
    shape_values[1][1, 2, 3] = np.float32(0.0)
    shape_values[2][2, 2, 2] = np.float32(0.0)
    channels.extend(memberships)
    channels.extend(
        [
            np.asarray(0.25 + x / 16 + y / 32 + z / 64, dtype=np.float32),
            np.asarray(0.125 + x / 32 + y / 64 + z / 128, dtype=np.float32),
            np.asarray(x / 16 + y / 32 + z / 64, dtype=np.float32),
        ]
    )
    channels.extend(class_data)
    channels.extend(
        [
            np.asarray(0.625 + x / 32 + y / 64 + z / 128, dtype=np.float32),
            np.asarray(0.5 + x / 32 + y / 64 + z / 128, dtype=np.float32),
            np.asarray(0.375 + x / 32 + y / 64 + z / 128, dtype=np.float32),
        ]
    )
    channels.extend(shape_values)
    channels.extend(
        [
            np.where(x == 1, 0, np.where(x == 2, 3, np.where(x == 3, -1, -1))),
            np.where(y == 1, 0, np.where(y == 2, 4, -1)),
            np.where(z == 2, 0, np.where(z == 3, 2, -1)),
        ]
    )
    return np.stack(channels, axis=-1).astype(np.float32)


def main() -> None:
    data = make_fixture()
    affine = np.diag([1.2, 1.4, 1.8, 1.0]).astype(np.float64)
    image = nib.Nifti1Image(data, affine)
    image.header.set_xyzt_units("mm")
    nib.save(image, HERE / "inputs_ras.nii.gz")


if __name__ == "__main__":
    main()
