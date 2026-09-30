from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
N = 33
NC = 2
K = 3
PAD = 1
BOUNDS = ((5, 36), (7, 38), (9, 40))
TRANSFORM = np.array([0.008, -0.012, 0.006, 0.08, -0.06, 0.05], dtype=np.float32)


def image_value(c: int, x: int, y: int, z: int) -> np.float32:
    if any(i < PAD or i >= N - PAD for i in (x, y, z)):
        return np.float32(0.0)
    fold = (x + 2 * y + 3 * z + c) % 7 - 3
    return np.float32((5 * x - 3 * y + 7 * z + 9 * c + fold) / 512.0)


images = np.zeros((N, N, N, NC), dtype=np.float32)
for x in range(N):
    for y in range(N):
        for z in range(N):
            for c in range(NC):
                images[x, y, z, c] = image_value(c, x, y, z)
image_path = ROOT / "padded_input_ras.nii.gz"
affine = np.diag([1.0, 1.0, 1.0, 1.0])
nib.save(nib.Nifti1Image(images, affine), image_path)

center = (N // 2, N // 2, N // 2)
centroids = np.array(
    [
        [image_value(0, *center), image_value(0, *center), -0.25],
        [image_value(1, *center), 0.35, -0.15],
    ],
    dtype=np.float32,
)

shape_dims = (35, 37, 36)
shape_files = []
for k in range(K):
    shape = np.empty(shape_dims, dtype=np.float32)
    for x in range(shape_dims[0]):
        for y in range(shape_dims[1]):
            for z in range(shape_dims[2]):
                fold = (3 * x + 5 * y + 7 * z + 11 * x * y + 13 * y * z + k) % 11 - 5
                value = np.float32(0.72 - 0.12 * k)
                value = np.float32(value + np.float32(0.0015 * x))
                value = np.float32(value + np.float32(0.0011 * y))
                value = np.float32(value + np.float32(0.0009 * z))
                shape[x, y, z] = np.float32(value + np.float32(fold * 0.0002))
    name = f"shape_{k}_ras.nii.gz"
    nib.save(nib.Nifti1Image(shape, affine), ROOT / name)
    shape_files.append(name)

metadata = {
    "input_nifti": image_path.name,
    "shape_niftis": shape_files,
    "input_shape_xyz_channels": list(images.shape),
    "shape_dimensions_xyz": list(shape_dims),
    "input_axis_codes": list(nib.aff2axcodes(affine)),
    "shape_axis_codes": list(nib.aff2axcodes(affine)),
    "centroids_channels_classes": centroids.tolist(),
    "crop_bounds_source_coordinates": [list(pair) for pair in BOUNDS],
    "pad_border": PAD,
    "preprocessing_transform": TRANSFORM.tolist(),
    "image_resolutions_mm": [1.0, 1.0, 1.0],
    "template_resolutions_mm": [1.0, 1.0, 1.0],
    "init_align_iterations": 50,
    "pyramid_levels": 4,
}
(ROOT / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
