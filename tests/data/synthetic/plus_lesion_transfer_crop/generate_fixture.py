from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (17, 17, 17)
AFFINE = np.array(
    [[0.8, 0.0, 0.0, -6.4], [0.0, 0.8, 0.0, -6.4], [0.0, 0.0, 1.2, -9.6], [0, 0, 0, 1]],
    dtype=np.float64,
)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    grid = np.indices(SHAPE, dtype=np.float32)
    x, y, z = grid[0], grid[1], grid[2]
    # An asymmetric folded tissue shell reaches the left image face so the
    # selected five-voxel crop border creates observable zero padding.
    center_x = x - 2.0
    center_y = y - 8.0
    center_z = z - 8.0
    radius = np.sqrt(center_x**2 + center_y**2 + center_z**2)
    theta = np.arctan2(center_y, center_x)
    fold = 0.7 * np.sin(3 * theta) + 0.25 * np.cos(2 * z)
    folded = radius + fold
    brain = folded < 7.0
    gm = (folded >= 3.4) & (folded < 6.0)
    wm = folded < 3.4

    sector = ((np.floor((theta + np.pi) * 2.0).astype(np.int16) % 5) + 200).astype(np.int16)
    seg = np.zeros(SHAPE, dtype=np.int16)
    seg[brain] = sector[brain]
    atlas_gm = np.where(gm, 175.0 + 0.25 * x + 0.1 * y, 0.0).astype(np.float32)
    atlas_wm = np.where(wm, 220.0, 0.0).astype(np.float32)
    atlas_mask = brain.astype(np.uint8) * 255
    brain_mask = (folded < 7.5).astype(np.uint8) * 255

    # An all-zero reference carries only the source grid geometry for the
    # selected P39 output arrays.
    arrays = {
        "lesion_output_grid.nii.gz": np.zeros(SHAPE, dtype=np.uint8),
        "segmentation_reference.nii.gz": seg,
        "multi_atlas_gm.nii.gz": atlas_gm,
        "multi_atlas_wm.nii.gz": atlas_wm,
        "multi_atlas_mask.nii.gz": atlas_mask,
        "brain_mask.nii.gz": brain_mask,
    }
    for name, values in arrays.items():
        image = nib.Nifti1Image(values, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, ROOT / name)


if __name__ == "__main__":
    main()
