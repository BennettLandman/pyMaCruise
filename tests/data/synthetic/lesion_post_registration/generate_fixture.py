from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SOURCE = ROOT.parent / "lesion_statistics_membership_bridge" / "bridge_inputs_ras.nii.gz"


def main() -> None:
    source = nib.load(str(SOURCE))
    base = np.asarray(source.dataobj, dtype=np.float32)
    shape = base.shape[:3]
    segmentation = np.zeros(shape, dtype=np.float32)
    # An asymmetric, folded class-1 ribbon with a class-2 interruption.
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                if (x + 2 * y + z) % 7 in (0, 1, 2) and x >= 2:
                    segmentation[x, y, z] = 1
                if x == 3 and y in (1, 2) and z == 2:
                    segmentation[x, y, z] = 2
    skeleton = np.full(shape, -1, dtype=np.float32)
    skeleton[2, 2, 1] = 1
    ordering = np.empty(shape, dtype=np.float32)
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                ordering[x, y, z] = -((13 * x + 7 * y + 3 * z) % 19) / 10.0
    frames = np.concatenate(
        (base, segmentation[..., None], skeleton[..., None], ordering[..., None]),
        axis=3,
    )
    output = nib.Nifti1Image(frames.astype(np.float32), source.affine, source.header)
    output.set_data_dtype(np.float32)
    nib.save(output, str(ROOT / "branch_inputs_ras.nii.gz"))


if __name__ == "__main__":
    main()
