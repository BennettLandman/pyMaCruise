#!/usr/bin/env python3
"""Create the fixed RAS state consumed by the P18czi branch adapter."""

from __future__ import annotations

import json
import struct
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/toads_post_registration_branch"
PREVIOUS_REFERENCE = ROOT / "tests/data/synthetic/toads_default_atlas_dispatch/java_reference.json"


def _save(name: str, data: np.ndarray, affine: np.ndarray) -> None:
    nib.save(nib.Nifti1Image(data, affine), str(FIXTURE / name))


def main() -> None:
    FIXTURE.mkdir(parents=True, exist_ok=True)
    shape = (9, 8, 7)
    affine = np.asarray(
        [[0.9, 0, 0, 12.5], [0, 1.2, 0, -31.0], [0, 0, 1.45, 8.25], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    reference = json.loads(PREVIOUS_REFERENCE.read_text())
    source_case = reference["cases"][5]  # no field/PSF, smoothed prior-bearing Java default sweep
    membership = np.asarray(source_case["membership_bits"], dtype=np.uint32).view(np.float32).reshape(*shape, 11)
    objective = struct.unpack("!f", struct.pack("!I", source_case["objective_bits"]))[0]

    x, y, z = np.indices(shape)
    # The asymmetric folded class-1 ribbon has two curved arms and a short
    # lateral branch; class-2 islands keep the labels non-concentric.
    center = 2.2 + 0.32 * y + 0.18 * np.sin(0.8 * z)
    folded = (y >= 1) & (y <= 6) & (z >= 1) & (z <= 5) & (np.abs(x - center) < 1.15)
    branch = (x >= 5) & (x <= 7) & (y >= 3) & (y <= 5) & (z >= 4) & (z <= 6)
    segmentation = np.zeros(shape, dtype=np.int16)
    segmentation[folded] = 1
    segmentation[branch] = 2

    skeleton = segmentation.astype(np.int16, copy=True)
    candidate = (4, 3, 3)
    segmentation[candidate] = 0
    skeleton[candidate] = -1
    # Ensure the synthetic missing skeleton point is surrounded by class 0,
    # so the selected grow initialization offers one deterministic event.
    segmentation[candidate[0]-1:candidate[0]+2, candidate[1]-1:candidate[1]+2, candidate[2]-1:candidate[2]+2] = 0
    skeleton[candidate[0]-1:candidate[0]+2, candidate[1]-1:candidate[1]+2, candidate[2]-1:candidate[2]+2] = 0
    skeleton[candidate] = -1

    # The source output transform is a fixed already-registered handoff. It is
    # recorded as metadata because the branch must carry it through unchanged.
    transform = np.asarray([0.018, -0.011, 0.007, 0.35, -0.22, 0.16], dtype=np.float32)
    ordering = (-0.01 * (x * 7 + y * 5 + z * 3)).astype(np.float32)
    class_volumes = np.asarray([18.0 + 0.25 * i for i in range(11)], dtype=np.float32)
    boundary_distances = np.zeros((*shape, 11), dtype=np.int8)
    template_labels = np.asarray([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], dtype=np.int8)
    previous_classification = np.zeros(shape, dtype=np.int8)

    for name, array in (
        ("membership_maps.nii.gz", membership),
        ("segmentation.nii.gz", segmentation),
        ("skeleton.nii.gz", skeleton),
        ("ordering.nii.gz", ordering),
        ("class_volumes.nii.gz", class_volumes),
        ("boundary_distances.nii.gz", boundary_distances),
        ("previous_classification.nii.gz", previous_classification),
    ):
        _save(name, array, affine)

    (FIXTURE / "state.json").write_text(
        json.dumps(
            {
                "source": "fixed state from selected 11-class Java reference sweep; registration output is supplied, not recomputed",
                "membership_source_case": 5,
                "membership_objective_bits": int(np.asarray([objective], np.float32).view(np.uint32)[0]),
                "registration_transform_bits": transform.view(np.uint32).astype(np.uint64).tolist(),
                "template_labels": template_labels.astype(int).tolist(),
                "grow_candidate_xyz": list(candidate),
                "true_ras_affine": affine.tolist(),
                "axes": "XYZ; RAS+",
            },
            indent=2,
        )
        + "\n"
    )
    folded_image = nib.load(str(FIXTURE / "segmentation.nii.gz"))
    assert nib.aff2axcodes(folded_image.affine) == ("R", "A", "S")
    assert segmentation.shape == (9, 8, 7)
    assert np.unique(segmentation).size == 3
    assert np.any(segmentation != segmentation[::-1, :, :])


if __name__ == "__main__":
    main()
