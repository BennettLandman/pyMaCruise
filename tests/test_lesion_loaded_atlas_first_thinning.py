from __future__ import annotations

from pathlib import Path

import numpy as np

from pymacruise.lesion_initial_atlas_alignment import (
    initialize_lesion_atlas_alignment_from_descriptor,
)
from pymacruise.lesion_loaded_atlas_first_thinning import (
    run_loaded_lesion_atlas_first_thinning,
)


ATLAS = (
    Path(__file__).parents[1]
    / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
)


def test_bundled_atlas_memberships_reach_first_topology_thinning_pass() -> None:
    """Exercise the loaded-atlas -> membership -> first thinning seam."""
    from pymacruise.lesion_atlas_priors import (
        load_lesion_atlas_priors,
        load_lesion_atlas_shape_priors,
    )

    priors = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(str(ATLAS))
    center = tuple(size // 2 for size in next(iter(shapes.shape_maps.values())).shape)
    start = tuple(value - 5 for value in center)
    slices = tuple(slice(value, value + 11) for value in start)
    atlas_patch = np.stack(
        [shapes.shape_maps[item.name][slices] for item in priors.structures], axis=-1
    ).astype(np.float32)
    dominant = np.argmax(atlas_patch, axis=-1)
    channels_last = priors.intensity_priors[:, dominant].transpose(1, 2, 3, 0)
    images = channels_last.transpose(3, 0, 1, 2).copy()
    aligned = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        images,
        np.zeros(6, dtype=np.float32),
        crop_bounds=((0, 9), (0, 9), (0, 9)),
        pad_border=1,
        image_resolutions_mm=(0.9375, 0.9375, 0.9375),
        template_resolutions_mm=(0.9375, 0.9375, 0.9375),
        init_align_iterations=0,
        pyramid_levels=2,
    )

    membership = aligned.alignment.memberships.memberships
    assert membership.shape == (11, 11, 11, 13)
    segmentation = np.argmax(membership, axis=-1).astype(np.int16)
    skeleton = segmentation.copy()
    available = np.ones((*segmentation.shape, 13), dtype=np.bool_)
    ordering = np.zeros(segmentation.shape, dtype=np.float32)
    auxiliary = np.zeros(segmentation.shape, dtype=np.float32)
    shape = np.stack(
        [aligned.shape_priors.shape_maps[item.name][slices] for item in priors.structures],
        axis=-1,
    ).astype(np.float32)
    object_types = tuple(item.topology for item in priors.structures)
    lesion_class = next(
        index for index, item in enumerate(priors.structures) if item.topology == "lesion"
    )

    result = run_loaded_lesion_atlas_first_thinning(
        aligned,
        segmentation,
        shape,
        auxiliary,
        auxiliary,
        skeleton,
        available,
        ordering,
        object_types=object_types,
        lesion_class=lesion_class,
        volume_ratio=0.0,
        limit=10.0,
    )

    assert result.alignment is aligned
    assert result.stage_order == (
        *aligned.alignment.call_order,
        "setGMflag",
        "initThinningLabels",
        "propagateThinning",
    )
    assert result.thinning.initialization.inserted_items
    assert result.thinning.run.events
    assert result.thinning.run.segmentation.shape == segmentation.shape
    assert np.isfinite(result.thinning.run.ordering).all()


def test_loaded_atlas_first_thinning_rejects_mismatched_topology_grid() -> None:
    from pymacruise.lesion_atlas_priors import (
        load_lesion_atlas_priors,
        load_lesion_atlas_shape_priors,
    )

    priors = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(str(ATLAS))
    center = tuple(size // 2 for size in next(iter(shapes.shape_maps.values())).shape)
    start = tuple(value - 5 for value in center)
    slices = tuple(slice(value, value + 11) for value in start)
    patch = np.stack(
        [shapes.shape_maps[item.name][slices] for item in priors.structures], axis=-1
    ).astype(np.float32)
    classes = np.argmax(patch, axis=-1)
    channels_last = priors.intensity_priors[:, classes].transpose(1, 2, 3, 0)
    images = channels_last.transpose(3, 0, 1, 2).copy()
    aligned = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS), images, np.zeros(6, dtype=np.float32),
        crop_bounds=((0, 9), (0, 9), (0, 9)), pad_border=1,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
        init_align_iterations=0, pyramid_levels=2,
    )
    labels = np.zeros((11, 11, 11), dtype=np.int16)
    maps = np.zeros((11, 11, 11), dtype=np.float32)
    with np.testing.assert_raises_regex(ValueError, "loaded atlas membership grid"):
        run_loaded_lesion_atlas_first_thinning(
            aligned, labels[:10], np.zeros((11, 11, 11, 13), dtype=np.float32),
            maps, maps, labels, np.ones((11, 11, 11, 13), dtype=np.bool_), maps,
            object_types=tuple(item.topology for item in priors.structures),
            lesion_class=10, volume_ratio=0.0, limit=0.0,
        )
