"""Build the P39zzzzs-derived folded RAS registration handoff fixture."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
BRIDGE = ROOT.parent / "lesion_statistics_membership_bridge" / "java_reference.json"
DIMS = (6, 5, 4)


def _shapes() -> np.ndarray:
    x, y, z = np.indices(DIMS)
    fold = (7 * x + 5 * y + 3 * z + 2 * x * y + 4 * y * z) % 11
    maps = np.zeros((*DIMS, 4), dtype=np.float32)
    maps[..., 0] = np.float32(0.9)
    maps[..., 1] = np.where(fold < 6, np.float32(0.9), np.float32(0.02))
    maps[0, :, :, 2] = np.float32(0.95)
    maps[-1, :, :, 2] = np.float32(0.8)
    maps[:, 0, :, 2] = np.float32(0.7)
    maps[:, -1, :, 2] = np.float32(0.6)
    # Class 3 is intentionally empty to pin source bbox sentinels.
    return maps


def _membership_maps() -> dict[str, np.ndarray]:
    reference = json.loads(BRIDGE.read_text())
    output: dict[str, np.ndarray] = {}
    for name in ("no_field", "with_field"):
        bits = np.asarray(reference[name]["class_memberships_bits"], dtype=np.uint32)
        output[name] = bits.view(np.float32).reshape((*DIMS, 4)).copy()
    return output


def make_fixture() -> None:
    memberships = _membership_maps()
    shapes = _shapes()
    frames = np.concatenate((memberships["no_field"], memberships["with_field"], shapes), axis=-1)
    affine = np.diag([0.9, 1.2, 1.5, 1.0]).astype(np.float64)
    affine[:3, 3] = (-8.1, 5.4, 12.0)
    image = nib.Nifti1Image(frames, affine)
    image.set_data_dtype(np.float32)
    nib.save(image, ROOT / "registration_handoff_ras.nii.gz")
    for name, values in memberships.items():
        flat = values.reshape(-1)
        (ROOT / f"{name}_memberships.txt").write_text(
            "\n".join(format(float(value), ".9g") for value in flat) + "\n"
        )


if __name__ == "__main__":
    make_fixture()
