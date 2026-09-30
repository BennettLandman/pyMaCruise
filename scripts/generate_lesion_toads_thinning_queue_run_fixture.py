"""Generate the synthetic RAS volume used by the LesionTOADS queue runner."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_thinning_queue_run"


def main() -> None:
    # Frames: segmentation, skeleton, class-0..2 memberships, lesion and
    # black-hole memberships, class-0..2 shape priors, class availability,
    # and ordering. All eight queued voxels are interior to this 5^3 image.
    packed = np.zeros((5, 5, 5, 14), dtype=np.float32)
    packed[..., 2] = 0.5
    packed[..., 3] = 0.25
    packed[..., 4] = 0.125
    for x in range(1, 4):
        for y in range(1, 4):
            for z in range(1, 4):
                if (x, y, z) not in {(2, 2, 2), (1, 1, 1), (1, 3, 1)}:
                    packed[x, y, z, 10] = 1.0
    packed[3, 1, 3, 1] = 2.0  # non-segmentation skeleton sentinel for a volume-gated item
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "queue_run_inputs.nii.gz"))


if __name__ == "__main__":
    main()
