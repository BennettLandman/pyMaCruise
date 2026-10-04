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
from .lesion_initial_alignment_bounds import (
    LesionInitialAlignmentBounds,
    compute_lesion_initial_alignment_bounds,
)
from .lesion_initial_atlas_alignment import LesionLoadedAtlasInitialAlignment
from .lesion_initial_segmentation_state import (
    LesionInitialSegmentationState,
    run_loaded_template_first_thinning,
)
from .lesion_registered_shapes import transform_registered_lesion_shape_priors
from .lesion_transformed_template import transform_lesion_toads_template
from .lesion_loaded_atlas_first_thinning import LesionLoadedAtlasFirstThinning


@dataclass(frozen=True)
class LesionRegisteredInitialThinning:
    """Derived RAS atlas maps, constructor state, and initial thinning."""

    alignment: LesionLoadedAtlasInitialAlignment
    topology: LesionAtlasTopology
    registration_bounds: LesionInitialAlignmentBounds
    transformed_template: NDArray[np.int8]
    transformed_shapes: NDArray[np.float32]
    constructor_state: LesionInitialSegmentationState
    first_thinning: LesionLoadedAtlasFirstThinning
    stage_order: tuple[str, ...]
    excluded_wrapper_stages: tuple[str, ...]


def run_registered_lesion_initial_thinning(
    alignment: LesionLoadedAtlasInitialAlignment,
    descriptor_path: str | Path,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
) -> LesionRegisteredInitialThinning:
    """Resample the registered 13-class atlas and execute first thinning.

    ``alignment`` must come from the selected descriptor-based initial
    membership and rigid-registration API. The topology is reloaded as RAS
    NIfTI; class shape maps stay in the same descriptor order. Constructor
    lesion and black-hole memberships start at the selected Java source's
    zero values.
    """
    if not isinstance(alignment, LesionLoadedAtlasInitialAlignment):
        raise TypeError("alignment must be a descriptor-loaded LesionTOADS alignment")
    structures = alignment.priors.structures
    if len(structures) != 13 or alignment.shape_priors.structures != structures:
        raise ValueError("selected LesionTOADS initial alignment requires its ordered 13-class atlas")
    memberships = alignment.alignment.memberships.memberships
    if memberships.ndim != 4 or memberships.shape[-1] != len(structures):
        raise ValueError("initial memberships must match the loaded atlas class order")
    dimensions = tuple(int(value) for value in memberships.shape[:3])
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
    black_hole_membership = np.zeros(dimensions, dtype=np.float32)
    initial = run_loaded_template_first_thinning(
        alignment,
        transformed_template,
        transformed_shapes,
        black_hole_membership,
        object_types=object_types,
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
        constructor_state=initial.state,
        first_thinning=initial.first_thinning,
        stage_order=stages,
        excluded_wrapper_stages=(
            "preprocess.uncropImages",
            "preprocess.findCroppingBoundaries",
            "preprocess.cropImages",
            "preprocess.estimateInitialTransform",
        ),
    )
