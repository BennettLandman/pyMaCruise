from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
shape = (3, 2, 4)
classes = 4
memberships = np.full((*shape, classes), -1.0, dtype=np.float32)
best = np.empty(shape, dtype=np.int16)
for x in range(shape[0]):
    for y in range(shape[1]):
        for z in range(shape[2]):
            best[x, y, z] = (x + 2 * y + 3 * z) % classes

def set_linear(linear_x_fast: int, class_index: int, value: float) -> None:
    x = linear_x_fast % shape[0]
    y = (linear_x_fast // shape[0]) % shape[1]
    z = linear_x_fast // (shape[0] * shape[1])
    memberships[x, y, z, class_index] = np.float32(value)

for k in range(classes):
    set_linear(0, k, 0.0)
    set_linear(1, k, -1.0 - k)
    set_linear(2, k, np.nan)
set_linear(3, 2, np.inf)
set_linear(4, 1, 0.75)
set_linear(4, 2, 0.75)
set_linear(5, 3, 0.75)
set_linear(5, 2, 0.25)
set_linear(6, 0, np.nan)
set_linear(6, 1, 0.2)
set_linear(7, 1, np.inf)
set_linear(7, 3, np.inf)
set_linear(8, 0, 0.5)
set_linear(8, 1, 0.8)
set_linear(8, 2, 0.8)
set_linear(8, 3, 0.1)
best[0, 0, :] = [1, 2, 3, 0]

affine = np.array(
    [[1.25, 0, 0, 14.0], [0, 1.5, 0, -22.0], [0, 0, 2.0, 8.0], [0, 0, 0, 1]],
    dtype=np.float64,
)
nib.save(nib.Nifti1Image(memberships, affine), HERE / "memberships.nii.gz")
nib.save(nib.Nifti1Image(best, affine), HERE / "best_indices.nii.gz")
