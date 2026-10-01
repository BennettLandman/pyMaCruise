from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "tests"))

from pymacruise.lesion_repeated_iteration import run_lesion_toads_repeated_iterations
from test_lesion_continuing_iteration import TYPES, _compose, _load, _options


def main() -> None:
    source_frames, affine = _load()
    _, initial_branch, _ = _compose(source_frames, use_field=False)
    result = run_lesion_toads_repeated_iterations(
        initial_branch,
        source_frames[..., :2],
        source_frames[..., 12:16],
        source_frames[..., 12:16],
        shape_priors=tuple(source_frames[..., 12 + index] for index in range(4)),
        initial_transform=np.asarray([0.03, -0.02, 0.01, 0.1, -0.1, 0.05], dtype=np.float32),
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
        bridge_options=_options(source_frames, use_field=False),
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        classification=np.zeros(source_frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        object_types=TYPES,
        lesion_class=1,
        n_iter_max=2,
        max_distance=0.0,
        algorithm_mode="standard",
        relation_limit=1,
        init_align_iter=1,
    )
    if len(result.passes) != 2:
        raise RuntimeError("the pinned fixed-input run must supply two pass states")

    # Seven channels per pass: segmentation, four class memberships, lesion
    # membership, and black-hole membership. Java consumes this RAS NIfTI
    # directly and computes initThinningLabels volume bits independently.
    nx, ny, nz = source_frames.shape[:3]
    states = np.zeros((nx, ny, nz, 14), dtype=np.float32)
    for pass_index, pass_state in enumerate(result.passes):
        start = pass_index * 7
        states[..., start] = pass_state.starting_segmentation
        states[..., start + 1 : start + 5] = pass_state.starting_class_memberships
        states[..., start + 5] = pass_state.starting_lesion_membership
        states[..., start + 6] = pass_state.starting_black_hole_membership
    path = ROOT / "tests/data/synthetic/lesion_repeated_iteration/loop_states_ras.nii.gz"
    nib.save(nib.Nifti1Image(states, affine), str(path))


if __name__ == "__main__":
    main()
