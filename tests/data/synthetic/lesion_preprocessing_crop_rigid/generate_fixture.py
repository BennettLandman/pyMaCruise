"""Generate a folded, asymmetric RAS crop fixture with a nonidentity transform."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
NX, NY, NZ, NC = 11, 10, 9, 2
RES = (0.8, 1.2, 1.5)


def image_value(c: int, x: int, y: int, z: int) -> np.float32:
    if not (3 <= x <= 7 and 2 <= y <= 7 and 2 <= z <= 6) or (x + 2*y + 3*z + c) % 6 == 0:
        return np.float32(0.0)
    return np.float32(((x*11 + y*7 + z*3 + c*5) % 17 + 1) / 20.0)


images = np.zeros((NX, NY, NZ, NC), dtype=np.float32)
for x in range(NX):
    for y in range(NY):
        for z in range(NZ):
            for c in range(NC):
                images[x, y, z, c] = image_value(c, x, y, z)
topology = np.ones((NX, NY, NZ), dtype=np.uint8)
topology[1:5, 3:8, 1:6] = 2
topology[7:10, 1:6, 3:8] = 3
affine = np.diag([*RES, 1.0])
image_name, topology_name = "normalized_channels_ras.nii.gz", "topology_labels_ras.nii.gz"
nib.save(nib.Nifti1Image(images, affine), ROOT / image_name)
nib.save(nib.Nifti1Image(topology, affine), ROOT / topology_name)
metadata = {
    "image_nifti": image_name,
    "topology_nifti": topology_name,
    "image_shape_xyz_channels": list(images.shape),
    "topology_shape_xyz": list(topology.shape),
    "axis_codes": list(nib.aff2axcodes(affine)),
    "image_resolution_mm": list(RES),
    "template_resolution_mm": list(RES),
    "orientation_signs": [1.0, 1.0, 1.0],
    "initial_transform": [0.04, -0.03, 0.08, -1.1, -2.6, -1.5],
    "background_label": 1,
    "background_ratio": 0.0,
    "padding": 10,
}
(ROOT / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
