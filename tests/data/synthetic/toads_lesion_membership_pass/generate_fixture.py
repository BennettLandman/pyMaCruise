"""Generate the deterministic folded base-TOADS LESION RAS fixture."""

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parent / "inputs_ras.nii.gz"
SPATIAL = (5, 4, 5)


def fixture() -> np.ndarray:
    nx, ny, nz = SPATIAL
    data = np.zeros((*SPATIAL, 16), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                code = x + 2 * y + 3 * z
                for c in range(2):
                    data[x, y, z, c] = np.float32(((7*x + 5*y + 3*z + 4*c) % 17 - 8) / 9.0)
                    field_delta = np.float32(np.float32((x + 2*y + z + c) % 5) * np.float32(0.09))
                    data[x, y, z, 2+c] = np.float32(np.float32(0.82) + field_delta)
                    data[x, y, z, 4+c] = np.float32(((2*x + 3*y + 5*z + c) % 9) / 31.0)
                if (x, y, z) == (2, 1, 2):
                    data[x, y, z, 0:2] = (np.float32(0.125), np.float32(-0.25))
                elif (x, y, z) == (3, 1, 2):
                    data[x, y, z, 0:2] = (np.float32(0.25), np.float32(-0.125))
                for k in range(2):
                    data[x, y, z, 6+k] = np.float32(((5*x + 7*y + 11*z + 3*k) % 31) / 30.0)
                    shape_delta = np.float32(np.float32((x + 3*y + 2*z + 5*k) % 11) * np.float32(0.04))
                    data[x, y, z, 9+k] = np.float32(np.float32(0.1) + shape_delta)
                data[x, y, z, 8] = np.float32(((3*x + 5*y + 7*z) % 13) / 16.0)
                data[x, y, z, 11] = np.float32(1 if code % 3 == 1 else 0)
                slots = [-1, -1, -1, -1]
                relation = code % 4
                if relation == 0:
                    slots[0] = 0
                    slots[3] = 0  # direct takes precedence over this neighbor slot
                elif relation == 1:
                    slots[1] = 1
                    slots[2] = 1  # LESION direct takes precedence over this neighbor slot
                elif relation == 2:
                    slots[2] = 1
                for slot, value in enumerate(slots):
                    data[x, y, z, 12+slot] = np.float32(value)
    return data


def main() -> None:
    array = fixture()
    image = nib.Nifti1Image(array, np.asarray([[1.2, 0, 0, 11], [0, 1.5, 0, -4], [0, 0, 2.1, 8], [0, 0, 0, 1]], dtype=np.float64))
    image.set_data_dtype(np.float32)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, str(OUTPUT))


if __name__ == "__main__":
    main()
