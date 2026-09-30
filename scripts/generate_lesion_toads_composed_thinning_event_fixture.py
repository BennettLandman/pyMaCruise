"""Generate the synthetic whole-volume inputs for one LesionTOADS event."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_composed_thinning_event"
NEIGHBORHOOD = "212221122120200120112112212"
N = 5
C = 3


def main() -> None:
    # Feature frames: segmentation, skeleton, three class memberships,
    # lesion map, black-hole map, three transformed-shape maps, three
    # availability maps, and the current ordering map.
    packed = np.zeros((N, N, N, 14), dtype=np.float32)
    packed[..., 0] = 0
    packed[..., 1] = 0
    for class_index, value in enumerate((0.125, 0.25, 0.0625)):
        packed[..., 2 + class_index] = value
    packed[..., 5] = 0.5
    packed[..., 6] = 0.25
    for index, char in enumerate(NEIGHBORHOOD):
        x, y, z = index // 9 + 1, (index // 3) % 3 + 1, index % 3 + 1
        packed[x, y, z, 0] = int(char)
        packed[x, y, z, 1] = int(char)
    packed[2, 2, 2, 1] = 0
    faces = ((1, 2, 2), (3, 2, 2), (2, 1, 2), (2, 3, 2), (2, 2, 1), (2, 2, 3))
    for x, y, z in faces:
        packed[x, y, z, 10] = 1.0
    packed[2, 2, 2, 13] = 0.5
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "event_inputs.nii.gz"))


if __name__ == "__main__":
    main()
