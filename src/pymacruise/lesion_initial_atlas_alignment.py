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
from .lesion_atlas_priors import (
    LesionAtlasPriors,
    LesionAtlasShapePriors,
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
)


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


@dataclass(frozen=True)
class LesionLoadedAtlasInitialAlignment:
    """Initial membership and registration using the bundled 13-class atlas.

    The descriptor supplies class-order intensity means and the separate shape
    loader supplies its matching RAS priors. Keeping both loaded products in
    the result makes their shared order available to subsequent selected
    thinning/registration stages without inventing a remapping.
    """

    priors: LesionAtlasPriors
    shape_priors: LesionAtlasShapePriors
    alignment: LesionInitialAtlasAlignment


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


def initialize_lesion_atlas_alignment_from_descriptor(
    descriptor_path: str,
    padded_images: ArrayLike,
    preprocessing_transform: ArrayLike,
    *,
    crop_bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    pad_border: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    init_align_iterations: int = 50,
    pyramid_levels: int = 4,
) -> LesionLoadedAtlasInitialAlignment:
    """Load the selected T1/FLAIR 13-class atlas and initialize its alignment.

    ``padded_images`` is the already normalized, source-ordered T1_MPRAGE and
    FLAIR pair with shape ``(2, x, y, z)``. The atlas descriptor provides the
    class centroids, and all thirteen shape maps are loaded in that same class
    order before invoking the source-matched initial-membership and rigid
    registration adapters. Image normalization, crop construction, transform
    estimation, topology initialization, and thinning are intentionally not
    synthesized here; those are separate upstream/next stages.
    """
    priors = load_lesion_atlas_priors(descriptor_path, ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(descriptor_path)
    if priors.structures != shapes.structures:
        raise ValueError("loaded intensity and shape atlases do not share class order")
    if len(priors.structures) != 13:
        raise ValueError("selected LesionTOADS alignment requires the 13-class atlas")
    data = np.asarray(padded_images, dtype=np.float32)
    if data.ndim != 4 or data.shape[0] != 2:
        raise ValueError("padded_images must have shape (2, x, y, z) in T1_MPRAGE/FLAIR order")
    alignment = initialize_lesion_atlas_alignment(
        data,
        priors.intensity_priors,
        tuple(shapes.shape_maps[item.name] for item in priors.structures),
        preprocessing_transform,
        crop_bounds=crop_bounds,
        pad_border=pad_border,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        init_align_iterations=init_align_iterations,
        pyramid_levels=pyramid_levels,
    )
    return LesionLoadedAtlasInitialAlignment(priors, shapes, alignment)
