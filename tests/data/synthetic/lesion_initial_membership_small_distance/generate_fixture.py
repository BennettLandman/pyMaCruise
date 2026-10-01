from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
PAD = 1
SHAPE = (4, 4, 4)
BOUNDS = ((5, 7), (8, 10), (2, 4))


def voxel(c: int, x: int, y: int, z: int) -> np.float32:
    if any(index < PAD or index >= SHAPE[axis] - PAD for axis, index in enumerate((x, y, z))):
        return np.float32(0.0)
    pattern = (x + 2 * y + 3 * z) % 3
    if pattern == 1:
        return np.float32(1.0e-12 if c in (0, 1) else 0.0)
    if pattern == 2 and c == 0:
        return np.float32(5.0e-11)
    return np.float32(0.0)


values = np.empty((*SHAPE, 2), dtype=np.float32)
for x in range(SHAPE[0]):
    for y in range(SHAPE[1]):
        for z in range(SHAPE[2]):
            for channel in range(2):
                values[x, y, z, channel] = voxel(channel, x, y, z)
affine = np.diag([0.8, 1.1, 1.5, 1.0])
path = ROOT / "small_distance_ras.nii.gz"
nib.save(nib.Nifti1Image(values, affine), path)
centroids = np.asarray([[0.0, 1.0e-12, 1.0e-9], [0.0, 1.0e-12, 0.0]], dtype=np.float32)
metadata = {
    "nifti": path.name,
    "shape_xyz_channels": list(values.shape),
    "crop_bounds_source_coordinates": [list(pair) for pair in BOUNDS],
    "pad_border": PAD,
    "centroids_channels_classes": centroids.tolist(),
    "centroid_bits": [int(value.view(np.uint32)) for value in centroids.ravel()],
    "source_zero": 1.0e-20,
}
(ROOT / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
