"""Source-ordered post-pyramid shape bounds for initial LesionTOADS alignment.

``LesionToadAtlas.registerShapesPyramid()`` finishes by refreshing positive
support bounds for each transformed atlas shape. This adapter composes that
handoff from an existing initial-alignment result; it does not repeat the
membership initialization or registration optimizer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_initial_atlas_alignment import LesionInitialAtlasAlignment
from .lesion_registration_handoff import _bounds_for_transform
from .toads_shapes import ToadShapeBounds


@dataclass(frozen=True)
class LesionInitialAlignmentBounds:
    """Final transform and per-shape bounds after the initial pyramid call."""

    transform: np.ndarray
    shape_bounds: tuple[ToadShapeBounds, ...]
    call_order: tuple[str, ...]


def compute_lesion_initial_alignment_bounds(
    alignment: LesionInitialAtlasAlignment,
    shape_priors: Sequence[ArrayLike],
    *,
    image_dimensions_xyz: tuple[int, int, int] | None = None,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
) -> LesionInitialAlignmentBounds:
    """Reproduce the selected post-pyramid transformed-shape bounds refresh.

    The result consumes the final transform from
    ``initialize_lesion_atlas_alignment``. By default it uses that result's
    subject grid; ``image_dimensions_xyz`` overrides the grid after a
    post-registration recrop. Bounds use the selected positive
    interpolated-support scan: minima start at the image dimensions, maxima
    start at zero, and upper coordinates are inclusive source maxima.
    """

    if not isinstance(alignment, LesionInitialAtlasAlignment):
        raise TypeError("alignment must be a LesionInitialAtlasAlignment result")
    transform = np.asarray(alignment.registration.transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("alignment must contain a finite six-value rigid transform")
    if len(shape_priors) != alignment.memberships.memberships.shape[-1]:
        raise ValueError("one shape prior is required for each initialized class")
    dimensions = (
        tuple(int(size) for size in alignment.memberships.memberships.shape[:3])
        if image_dimensions_xyz is None
        else tuple(int(size) for size in image_dimensions_xyz)
    )
    if len(dimensions) != 3 or any(size < 2 for size in dimensions):
        raise ValueError("image_dimensions_xyz must contain three dimensions >= 2")
    bounds = tuple(
        _bounds_for_transform(
            shape,
            image_dimensions_xyz=dimensions,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=transform,
        )
        for shape in shape_priors
    )
    return LesionInitialAlignmentBounds(
        transform=transform.copy(),
        shape_bounds=bounds,
        call_order=("atlas.registerShapesPyramid", "atlas.computeTransformedShapeBoundingBox"),
    )
