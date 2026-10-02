"""Generate the sparse integer MaSeg labels used by the base workflow run."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE_T1 = ROOT / "toads_registered_first_iteration_cycle" / "t1.nii.gz"
OUTPUT = Path(__file__).resolve().parent / "maseg.nii.gz"
LABELS = (35, 38, 39, 40, 41, 71, 72, 73)
COORDINATES = (
    (3, 3, 3), (3, 3, 9), (3, 9, 3), (3, 9, 9),
    (9, 3, 3), (9, 3, 9), (9, 9, 3), (9, 9, 9),
)


def main() -> None:
    source = nib.load(str(SOURCE_T1))
    labels = np.zeros(source.shape, dtype=np.uint8)
    for coordinate, label in zip(COORDINATES, LABELS, strict=True):
        labels[coordinate] = label
    output = nib.Nifti1Image(labels, source.affine, header=source.header.copy())
    output.set_data_dtype(np.uint8)
    nib.save(output, str(OUTPUT))
    print(OUTPUT)


if __name__ == "__main__":
    main()
