"""Continue bundled-atlas LesionTOADS alignment through first thinning.

This bounded composition consumes the already initialized descriptor
membership/rigid-registration result, resamples its topology and shape maps in
RAS voxel order, creates selected constructor state, and starts the first
topology-preserving thinning pass. It does not perform input normalization,
cropping, alignment-transform seeding, or later LesionTOADS iterations.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .lesion_atlas_priors import LesionAtlasTopology, load_lesion_atlas_topology
from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_initial_alignment_bounds import (
    LesionInitialAlignmentBounds,
    compute_lesion_initial_alignment_bounds,
)
from .lesion_initial_atlas_alignment import LesionLoadedAtlasInitialAlignment
from .lesion_initial_segmentation_state import (
    LesionInitialSegmentationState,
    initialize_lesion_segmentation_state,
)
from .lesion_registered_shapes import transform_registered_lesion_shape_priors
from .lesion_transformed_template import transform_lesion_toads_template
from .lesion_standard13_volume import Standard13VolumeResult


@dataclass(frozen=True)
class LesionRegisteredInitialThinning:
    """Derived RAS atlas maps, constructor state, and initial thinning."""

    alignment: LesionLoadedAtlasInitialAlignment
    topology: LesionAtlasTopology
    registration_bounds: LesionInitialAlignmentBounds
    transformed_template: NDArray[np.int8]
    transformed_shapes: NDArray[np.float32]
    constructor_state: LesionInitialSegmentationState
    initial_membership: Standard13VolumeResult
    first_thinning: LesionToadsThinningPass
    stage_order: tuple[str, ...]
    excluded_wrapper_stages: tuple[str, ...]


def run_registered_lesion_initial_thinning(
    alignment: LesionLoadedAtlasInitialAlignment,
    descriptor_path: str | Path,
    *,
    initial_membership: Standard13VolumeResult,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
) -> LesionRegisteredInitialThinning:
    """Resample the registered 13-class atlas and execute first thinning.

    ``alignment`` contains the pre-registration ``preprocess.initialMemberships``
    consumed by atlas alignment. ``initial_membership`` is a separate result
    from the constructor's later ``computeMemberships`` call: Java discards the
    former after registration and does not reuse it as segmentation weights.
    This adapter consumes (but does not calculate) that post-constructor
    result. The topology is reloaded as RAS NIfTI; class shape maps stay in the
    descriptor order.
    """
    if not isinstance(alignment, LesionLoadedAtlasInitialAlignment):
        raise TypeError("alignment must be a descriptor-loaded LesionTOADS alignment")
    if not isinstance(initial_membership, Standard13VolumeResult):
        raise TypeError("initial_membership must be a selected 13-class computeMemberships result")
    structures = alignment.priors.structures
    if len(structures) != 13 or alignment.shape_priors.structures != structures:
        raise ValueError("selected LesionTOADS initial alignment requires its ordered 13-class atlas")
    memberships = alignment.alignment.memberships.memberships
    if memberships.ndim != 4 or memberships.shape[-1] != len(structures):
        raise ValueError("initial memberships must match the loaded atlas class order")
    dimensions = tuple(int(value) for value in memberships.shape[:3])
    if initial_membership.class_memberships.shape != (*dimensions, 13):
        raise ValueError("post-constructor memberships must match the transformed subject grid")
    atlas_shapes = tuple(
        alignment.shape_priors.shape_maps[item.name] for item in structures
    )
    bounds = compute_lesion_initial_alignment_bounds(
        alignment.alignment,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
    )
    topology = load_lesion_atlas_topology(descriptor_path)
    template_values = np.asarray(topology.image.dataobj)
    transformed_template = transform_lesion_toads_template(
        template_values,
        bounds.transform,
        dimensions,
        image_resolution_mm=image_resolutions_mm,
        template_resolution_mm=template_resolutions_mm,
        background_label=structures[0].label,
    )
    transformed_shapes = transform_registered_lesion_shape_priors(
        bounds,
        atlas_shapes,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
    )
    object_types = tuple(item.topology for item in structures)
    lesion_classes = [index for index, item in enumerate(structures) if item.topology.lower() == "lesion"]
    if len(lesion_classes) != 1:
        raise ValueError("selected LesionTOADS descriptor must declare exactly one lesion class")
    state = initialize_lesion_segmentation_state(
        transformed_template,
        tuple(item.label for item in structures),
    )
    if state.segmentation.shape != initial_membership.class_memberships.shape[:3]:
        raise ValueError("transformed template and post-constructor memberships must share a grid")
    thinning = propagate_lesion_toads_thinning(
        state.segmentation,
        initial_membership.class_memberships,
        object_types,
        initial_membership.lesion_membership,
        initial_membership.black_hole_membership,
        transformed_shapes,
        state.skeleton,
        state.available,
        state.ordering,
        lesion_class=lesion_classes[0],
        volume_ratio=volume_ratio,
        limit=limit,
        use_lesion_max=use_lesion_max,
    )
    stages = (
        *alignment.alignment.call_order,
        "atlas.computeTransformedShapeBoundingBox",
        "atlas.computeTransformedTemplate.afterInitialRegistration",
        "atlas.getFastTransformedShape",
        "MultipleToadSegmentation.constructor:initializeClassState",
        "LesionToadSegmentation.computeMemberships:GMflag=false (supplied result)",
        "setGMflag",
        "initThinningLabels",
        "propagateThinning",
    )
    return LesionRegisteredInitialThinning(
        alignment=alignment,
        topology=topology,
        registration_bounds=bounds,
        transformed_template=transformed_template,
        transformed_shapes=transformed_shapes,
        constructor_state=state,
        initial_membership=initial_membership,
        first_thinning=thinning,
        stage_order=stages,
        excluded_wrapper_stages=(
            "preprocess.uncropImages",
            "preprocess.findCroppingBoundaries",
            "preprocess.cropImages",
            "preprocess.estimateInitialTransform",
        ),
    )
