"""Selected LesionTOADS initial membership-to-atlas registration handoff."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_initial_memberships import (
    LesionInitialMemberships,
    lesion_initial_memberships,
)
from .lesion_shapes import LesionRegistrationPyramidResult, register_lesion_shapes_pyramid


_CALL_ORDER = (
    "preprocess.initialMemberships",
    'atlas.initShapeRegistration(memberships,"rigid",initAlignIter,levels,1.0f)',
    "atlas.setTransform(preprocess.getTransform())",
    'atlas.updateShapeRegistration(memberships,"rigid",initAlignIter,levels,1.0f)',
    "atlas.registerShapesPyramid",
)


@dataclass(frozen=True)
class LesionInitialAtlasAlignment:
    """The exact initializer map plus selected initial pyramid registration."""

    memberships: LesionInitialMemberships
    registration: LesionRegistrationPyramidResult
    preprocessing_transform: NDArray[np.float32]
    call_order: tuple[str, ...]
    initialization_iterations: int
    pyramid_levels: int
    registration_scale_factor: int
    coupling_factor: np.float32


def initialize_lesion_atlas_alignment(
    padded_images: ArrayLike,
    centroids: ArrayLike,
    shape_priors: Sequence[ArrayLike],
    preprocessing_transform: ArrayLike,
    *,
    crop_bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    pad_border: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    init_align_iterations: int = 50,
    pyramid_levels: int = 4,
) -> LesionInitialAtlasAlignment:
    """Compose the active wrapper's initial rigid alignment sequence.

    This joins selected `initialMemberships()` to `initShapeRegistration`,
    the preprocessing-seeded transform, same-mode `updateShapeRegistration`,
    and `registerShapesPyramid`. The two registration setup calls both select
    rigid mode, `subsample=2`, and coupling factor 1.0; the already validated
    P39o pyramid optimizer performs the actual scale loop. Input arrays follow
    RAS voxel-axis order. Cropping, preprocessing-transform calculation,
    atlas image loading, and later single-scale segmentation alignment are
    outside this composition.
    """
    if (
        isinstance(init_align_iterations, (bool, np.bool_))
        or not isinstance(init_align_iterations, (int, np.integer))
        or init_align_iterations < 0
    ):
        raise ValueError("init_align_iterations must be a nonnegative integer")
    if isinstance(pyramid_levels, (bool, np.bool_)) or not isinstance(
        pyramid_levels, (int, np.integer)
    ) or pyramid_levels < 1:
        raise ValueError("pyramid_levels must be a positive integer")
    transform = np.asarray(preprocessing_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("preprocessing_transform must contain six finite rigid values")
    transform = transform.copy()

    initialized = lesion_initial_memberships(
        padded_images,
        centroids,
        crop_bounds=crop_bounds,
        pad_border=pad_border,
    )
    registration = register_lesion_shapes_pyramid(
        shape_priors,
        initialized.memberships,
        transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        levels=int(pyramid_levels),
        subsample=2,
        max_iterations=int(init_align_iterations),
    )
    return LesionInitialAtlasAlignment(
        memberships=initialized,
        registration=registration,
        preprocessing_transform=transform,
        call_order=_CALL_ORDER,
        initialization_iterations=int(init_align_iterations),
        pyramid_levels=int(pyramid_levels),
        registration_scale_factor=2,
        coupling_factor=np.float32(1.0),
    )
