"""Generate deterministic RAS NIfTI inputs for the LesionTOADS queue loop."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_thinning_schedule"
CASES = (
    ("all_faces_ordinary", 0),
    ("lesion_max", 1),
    ("lesion_sum", 1),
    ("ratio_tie_limit_equal", 0),
    ("ratio_tie_limit_above", 0),
    ("volume_equal", 0),
    ("volume_below", 0),
    ("availability_gate", 0),
)
FACE = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))


def main() -> None:
    packed = np.zeros((5, 5, 5, len(CASES) * 6), dtype=np.float32)
    for case_index, (name, removed_class) in enumerate(CASES):
        base = case_index * 6
        for x in range(5):
            for y in range(5):
                for z in range(5):
                    for c in range(3):
                        packed[x, y, z, base + c] = (2 + ((x * 7 + y * 5 + z * 3 + c * 11) % 15)) / 32.0
                    packed[x, y, z, base + 3] = (1 + ((x * 3 + y * 5 + z * 7) % 9)) / 32.0
                    packed[x, y, z, base + 4] = (1 + ((x * 5 + y * 3 + z * 2) % 7)) / 64.0
        if name in {"ratio_tie_limit_equal", "ratio_tie_limit_above"}:
            for dx, dy, dz in FACE:
                x, y, z = (2 + dx, 2 + dy, 2 + dz)
                packed[x, y, z, base : base + 3] = 0.25
                packed[x, y, z, base + 3 : base + 5] = 0.0
        for dx, dy, dz in FACE:
            x, y, z = (2 + dx, 2 + dy, 2 + dz)
            packed[x, y, z, base + 5] = 1.0
        if name == "availability_gate":
            packed[1, 2, 2, base + 5] = 0.0
        if removed_class not in range(3):
            raise AssertionError("bad class index")
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "schedule_inputs.nii.gz"))


if __name__ == "__main__":
    main()
