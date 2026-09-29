"""Write the deterministic RAS label volume used by P39s."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
N = 9
AFFINE = np.diag([1.5, 2.0, 2.5, 1.0])


def make_labels() -> np.ndarray:
    labels = np.zeros((N, N, N), dtype=np.int16)
    for x in range(N):
        for y in range(N):
            for z in range(N):
                center = 2 + ((y + z) % 3)
                if 2 <= y <= 6 and 2 <= z <= 6 and abs(x - center) <= 1:
                    labels[x, y, z] = 1
                if 5 <= x <= 6 and 2 <= y <= 3 and 5 <= z <= 6:
                    labels[x, y, z] = 2
                if x == 1 and 2 <= y <= 6 and z == 7:
                    labels[x, y, z] = 3
    return labels


if __name__ == "__main__":
    image = nib.Nifti1Image(make_labels(), AFFINE)
    image.set_qform(AFFINE, code=1)
    image.set_sform(AFFINE, code=1)
    nib.save(image, ROOT / "folded_labels.nii.gz")
