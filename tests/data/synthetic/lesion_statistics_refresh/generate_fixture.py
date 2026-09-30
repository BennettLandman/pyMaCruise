from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (6, 5, 4)
CLASSES = 5
CHANNELS = 2


def make_frames() -> np.ndarray:
    frames = np.zeros((*SHAPE, 1 + CLASSES + 2 * CHANNELS + 3), dtype=np.float32)
    for x in range(SHAPE[0]):
        for y in range(SHAPE[1]):
            for z in range(SHAPE[2]):
                label = (x + 2 * y + z + ((x * y + 2 * z) % 3)) % 4
                frames[x, y, z, 0] = label
                for k in range(CLASSES):
                    frames[x, y, z, 1 + k] = np.float32(
                        np.float32(0.08)
                        + np.float32(np.float32((x * 3 + y * 5 + z * 7 + k * 2) % 9) * np.float32(0.055))
                    )
                for c in range(CHANNELS):
                    frames[x, y, z, 1 + CLASSES + c] = np.float32(
                        np.float32(0.25)
                        + np.float32(np.float32((x * 11 + y * 7 + z * 3 + c * 5) % 23) * np.float32(0.12))
                    )
                    frames[x, y, z, 1 + CLASSES + CHANNELS + c] = np.float32(
                        np.float32(0.75)
                        + np.float32(np.float32((x * 2 + y * 3 + z * 4 + c) % 5) * np.float32(0.11))
                    )
                frames[x, y, z, 1 + CLASSES + 2 * CHANNELS] = np.float32(
                    np.float32((x + 3 * y + z) % 5) * np.float32(0.13)
                )
                frames[x, y, z, 1 + CLASSES + 2 * CHANNELS + 1] = np.float32(
                    np.float32(0.04)
                    + np.float32(np.float32((x + 2 * y + z) % 9) * np.float32(0.07))
                )
                frames[x, y, z, 1 + CLASSES + 2 * CHANNELS + 2] = np.float32(
                    np.float32(0.03)
                    + np.float32(np.float32((2 * x + y + 3 * z) % 9) * np.float32(0.075))
                )

    # Pin source-specific LESION and BLACKHOLE winner/tie branches.
    points: dict[int, list[tuple[int, int, int]]] = {k: [] for k in range(4)}
    for x in range(SHAPE[0]):
        for y in range(SHAPE[1]):
            for z in range(SHAPE[2]):
                points[int(frames[x, y, z, 0])].append((x, y, z))
    lesion = 1 + CLASSES + 2 * CHANNELS + 1
    black = lesion + 1
    for point, k, m, lm, bm in (
        (points[2][0], 2, 0.30, 0.80, 0.20),
        (points[2][1], 2, 0.25, 0.70, 0.70),
        (points[3][0], 3, 0.25, 0.20, 0.90),
    ):
        x, y, z = point
        frames[x, y, z, 1 + k] = np.float32(m)
        frames[x, y, z, lesion] = np.float32(lm)
        frames[x, y, z, black] = np.float32(bm)
    return frames


if __name__ == "__main__":
    data = make_frames()
    nii = nib.Nifti1Image(data, np.eye(4, dtype=np.float64))
    nii.set_qform(np.eye(4), code=1)
    nii.set_sform(np.eye(4), code=1)
    nib.save(nii, str(ROOT / "statistics_refresh_inputs.nii.gz"))
