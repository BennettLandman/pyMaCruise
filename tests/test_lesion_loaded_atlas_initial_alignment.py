from __future__ import annotations

from pathlib import Path

import numpy as np

from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
)
from pymacruise.lesion_initial_atlas_alignment import (
    initialize_lesion_atlas_alignment_from_descriptor,
)
from pymacruise.lesion_initial_memberships import lesion_initial_memberships


ATLAS = (
    Path(__file__).parents[1]
    / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
)


def test_bundled_13_class_atlas_and_synthetic_normalized_t1_flair_reach_initial_alignment():
    priors = load_lesion_atlas_priors(ATLAS, ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(ATLAS)
    assert priors.structures == shapes.structures
    assert len(priors.structures) == 13

    # Build a compact normalized two-channel phantom from the real atlas patch:
    # each voxel receives its dominant class's descriptor intensity. Two pad
    # layers are left for the selected active-crop membership operation.
    center = tuple(size // 2 for size in next(iter(shapes.shape_maps.values())).shape)
    starts = tuple(value - 5 for value in center)
    slices = tuple(slice(value, value + 11) for value in starts)
    patch = np.stack(
        [shapes.shape_maps[item.name][slices] for item in priors.structures], axis=-1
    ).astype(np.float32)
    classes = np.argmax(patch, axis=-1)
    channels_last = priors.intensity_priors[:, classes].transpose(1, 2, 3, 0)
    padded_images = np.transpose(channels_last, (3, 0, 1, 2)).copy()

    result = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        padded_images,
        np.zeros(6, dtype=np.float32),
        crop_bounds=((0, 9), (0, 9), (0, 9)),
        pad_border=1,
        image_resolutions_mm=(0.9375, 0.9375, 0.9375),
        template_resolutions_mm=(0.9375, 0.9375, 0.9375),
        init_align_iterations=0,
        pyramid_levels=2,
    )

    assert result.priors.structures == result.shape_priors.structures
    assert tuple(result.alignment.memberships.memberships.shape) == (11, 11, 11, 13)
    assert result.alignment.memberships.active_bounds == ((1, 10),) * 3
    assert tuple(result.alignment.registration.transform.shape) == (6,)
    assert result.alignment.call_order[0] == "preprocess.initialMemberships"
    assert result.alignment.call_order[-1] == "atlas.registerShapesPyramid"
    assert tuple(level.scale for level in result.alignment.registration.levels) == (2,)

    direct = lesion_initial_memberships(
        padded_images,
        result.priors.intensity_priors,
        crop_bounds=((0, 9), (0, 9), (0, 9)),
        pad_border=1,
    )
    np.testing.assert_array_equal(
        result.alignment.memberships.memberships.view(np.int32),
        direct.memberships.view(np.int32),
    )
