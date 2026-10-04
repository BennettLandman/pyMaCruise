"""Connect loaded-atlas LesionTOADS initialization to its first thinning pass.

This adapter consumes the initial memberships produced by the selected
``LesionToadPreprocess.initialMemberships`` path and carries those same maps
into ``initThinningLabels`` and ``propagateThinning``. Segmentation,
transformed shape, auxiliary lesion maps, skeleton, availability, and ordering
are explicit inputs because their full construction is not yet composed here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_initial_atlas_alignment import LesionLoadedAtlasInitialAlignment


@dataclass(frozen=True)
class LesionLoadedAtlasFirstThinning:
    """Initial atlas alignment and its first source-ordered thinning result."""

    alignment: LesionLoadedAtlasInitialAlignment
    thinning: LesionToadsThinningPass
    stage_order: tuple[str, ...]


def run_loaded_lesion_atlas_first_thinning(
    alignment: LesionLoadedAtlasInitialAlignment,
    segmentation: ArrayLike,
    transformed_shape: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    initial_skeleton: ArrayLike,
    initial_available: ArrayLike,
    initial_ordering: ArrayLike,
    *,
    object_types: Sequence[str],
    lesion_class: int,
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
) -> LesionLoadedAtlasFirstThinning:
    """Run one first-loop thinning pass from a descriptor-loaded atlas.

    The membership array is intentionally taken from ``alignment`` instead of
    accepting a second membership argument, so the loaded-atlas result cannot
    be accidentally bypassed at this workflow seam. Caller-supplied images
    must occupy its same RAS voxel grid. Earlier preprocessing, topology
    construction, and transformed-shape resampling remain upstream inputs.
    """
    if not isinstance(alignment, LesionLoadedAtlasInitialAlignment):
        raise TypeError("alignment must be a LesionLoadedAtlasInitialAlignment")
    memberships = alignment.alignment.memberships.memberships
    if memberships.ndim != 4 or memberships.shape[-1] != 13:
        raise ValueError("loaded LesionTOADS memberships must contain 13 classes")
    if len(object_types) != memberships.shape[-1]:
        raise ValueError("object_types must match the loaded 13-class atlas")
    spatial_shape = memberships.shape[:3]
    for name, value in (
        ("segmentation", segmentation),
        ("transformed_shape", transformed_shape),
        ("lesion_membership", lesion_membership),
        ("black_hole_membership", black_hole_membership),
        ("initial_skeleton", initial_skeleton),
        ("initial_available", initial_available),
        ("initial_ordering", initial_ordering),
    ):
        array = np.asarray(value)
        if array.shape[:3] != spatial_shape:
            raise ValueError(f"{name} must use the loaded atlas membership grid")

    thinning = propagate_lesion_toads_thinning(
        segmentation,
        memberships,
        object_types,
        lesion_membership,
        black_hole_membership,
        transformed_shape,
        initial_skeleton,
        initial_available,
        initial_ordering,
        lesion_class=lesion_class,
        volume_ratio=volume_ratio,
        limit=limit,
        use_lesion_max=use_lesion_max,
    )
    return LesionLoadedAtlasFirstThinning(
        alignment=alignment,
        thinning=thinning,
        stage_order=(
            *alignment.alignment.call_order,
            "setGMflag",
            "initThinningLabels",
            "propagateThinning",
        ),
    )
