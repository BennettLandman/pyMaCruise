from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
NX, NY, NZ, CHANNELS, PAD = 7, 8, 9, 2, 1
CROP_BOUNDS = ((2, 7), (3, 9), (4, 11))


def image(c: int, x: int, y: int, z: int) -> np.float32:
    if x < PAD or x >= NX - PAD or y < PAD or y >= NY - PAD or z < PAD or z >= NZ - PAD:
        return np.float32(0)
    fold = (x + 2 * y + 3 * z + c) % 5 - 2
    return np.float32((7 * x - 5 * y + 3 * z + 11 * c + fold) / 64.0)


data = np.zeros((NX, NY, NZ, CHANNELS), dtype=np.float32)
for x in range(NX):
    for y in range(NY):
        for z in range(NZ):
            for c in range(CHANNELS):
                data[x, y, z, c] = image(c, x, y, z)
affine = np.diag([1.0, 1.0, 1.0, 1.0])
image_file = ROOT / "padded_multichannel_ras.nii.gz"
nib.save(nib.Nifti1Image(data, affine), image_file)

centroids = np.array(
    [
        [image(0, 1, 1, 1), image(0, 1, 1, 1), -0.25, 0.125],
        [image(1, 1, 1, 1), 0.375, 0.4, -0.375],
    ],
    dtype=np.float32,
)
metadata = {
    "nifti": image_file.name,
    "shape_xyz_channels": list(data.shape),
    "axis_codes": list(nib.aff2axcodes(affine)),
    "crop_bounds_source_coordinates": [list(v) for v in CROP_BOUNDS],
    "pad_border": PAD,
    "centroids_channels_classes": centroids.tolist(),
    "centroid_bits": [int(value.view(np.uint32)) for value in centroids.ravel()],
}
(ROOT / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
