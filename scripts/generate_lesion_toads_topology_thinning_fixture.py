"""Generate a folded, mostly coherent synthetic LesionTOADS thinning case."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_topology_thinning"
N = 11


def _labels() -> np.ndarray:
    """Use the documented 31-bit LCG so Java and Python generate identical labels."""
    state = 1
    labels = np.empty((N, N, N), dtype=np.int16)
    for x in range(N):
        for y in range(N):
            for z in range(N):
                base = ((x // 3) + 2 * (y // 3) + (z // 3)) % 3
                state = (1103515245 * state + 12345) & 0x7FFFFFFF
                labels[x, y, z] = (state // 100) % 3 if state % 100 < 5 else base
    return labels


def main() -> None:
    # Frames: segmentation, skeleton, three memberships, lesion and black-hole
    # membership, three transformed shapes, three availability maps, ordering.
    packed = np.zeros((N, N, N, 14), dtype=np.float32)
    ROOT.mkdir(parents=True, exist_ok=True)
    labels = _labels()
    packed[..., 0] = labels
    x, y, z = np.indices(labels.shape)
    border = ((x == 0) | (x == N - 1) | (y == 0) | (y == N - 1)
              | (z == 0) | (z == N - 1))
    packed[..., 1] = np.where(border, 7, -1)
    for k in range(3):
        selected = labels == k
        packed[..., 2 + k][selected] = 0.49
        packed[..., 2 + (k + 1) % 3][selected] = 0.51
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "topology_thinning_inputs.nii.gz"))


if __name__ == "__main__":
    main()
