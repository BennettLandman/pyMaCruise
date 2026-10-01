"""Generate the asymmetric folded RAS inputs for the P40 pre-TGDM slice."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent


def main() -> None:
    size = 25
    grid = np.indices((size, size, size), dtype=np.float32)
    x, y, z = (grid[axis] - 12.0 for axis in range(3))
    radius = np.sqrt(x * x + y * y + z * z)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, np.sqrt(x * x + y * y))
    folds = 0.85 * np.sin(3.0 * theta + 0.4) * np.cos(2.0 * phi)
    folds += 0.35 * np.sin(5.0 * theta - 0.7) * np.cos(phi)
    folded_radius = radius + folds
    wm = folded_radius < 5.2
    gm = (folded_radius >= 5.2) & (folded_radius < 8.6)
    brain = folded_radius < 9.7

    sectors = np.floor((theta + np.pi) * 2.5).astype(np.int32) % 5
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    segmentation[brain] = (200 + sectors[brain]).astype(np.int16)
    segmentation[gm] = (210 + sectors[gm]).astype(np.int16)
    segmentation[gm & (x > 1.0) & (y < -1.0)] = 202

    atlas_gm = np.where(brain, 180.0, 0.0).astype(np.float32)
    atlas_wm = np.where(wm, 255.0, 0.0).astype(np.float32)
    toads_gm = np.where(gm, 200.0, 0.0).astype(np.float32)
    toads_wm = np.where(wm, 255.0, 0.0).astype(np.float32)
    brain_mask = brain.astype(np.uint8) * 255
    affine = np.diag((0.8, 0.8, 1.2, 1.0))

    arrays = {
        "segmentation.nii.gz": segmentation,
        "toads_wm_for_boundary.nii.gz": toads_wm.copy(),
        "multi_atlas_gm.nii.gz": atlas_gm,
        "multi_atlas_wm.nii.gz": atlas_wm,
        "multi_atlas_brain_mask.nii.gz": brain_mask,
        "toads_gm.nii.gz": toads_gm,
        "toads_wm.nii.gz": toads_wm,
        "toads_mask.nii.gz": brain_mask.copy(),
    }
    hashes = {}
    for name, values in arrays.items():
        path = ROOT / name
        nib.save(nib.Nifti1Image(values, affine), str(path))
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    provenance = {
        "kind": "project-generated synthetic folded-cortex stage fixture",
        "shape": [size, size, size],
        "affine_ras_mm": affine.tolist(),
        "sha256": hashes,
    }
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
