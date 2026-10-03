"""Continue a registered exact-membership update through its selected branch."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_alignment_startup import ToadsAtlasAlignmentRun
from .toads_atlas import toads_atlas_shapes_in_structure_order
from .toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from .toads_default_membership_registration_control import (
    ToadsDefaultMembershipRegistrationControl,
    continue_default_membership_through_registration_control,
)
from .toads_default_membership_registered_iteration import (
    ToadsDefaultMembershipRegisteredIteration,
    continue_default_membership_registered_iteration,
)
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)
from .toads_shape_volume import transform_toads_shape_priors


@dataclass(frozen=True)
class ToadsRegisteredMembershipBranchCycle:
    """One exact membership update plus registration and its wrapper branch."""

    membership_pass: ToadsDefaultAtlasMembershipPass
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    next_transformed_shape_priors: NDArray[np.float32] | None
    next_thinning: ToadsPostRegistrationNextThinning | None
    final_segmentation: NDArray
    stage_order: tuple[str, ...]


def continue_registered_membership_through_branch(
    alignment_run: ToadsAtlasAlignmentRun,
    membership_pass: ToadsDefaultAtlasMembershipPass,
    *,
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    boundary_distances: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    relation_limit: int,
    do_growing: bool,
    initial_transform: ArrayLike,
    initial_available: ArrayLike | None = None,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    record_growing_events: bool = True,
    record_thinning_items: bool = True,
) -> ToadsRegisteredMembershipBranchCycle:
    """Apply membership → registration → stop/grow-or-refresh → next thinning.

    This stateful adapter reuses the selected component implementations. It
    transforms shape priors from the *resulting* rigid pose before a
    continuing branch enters its next loop, avoiding caller-supplied shapes
    that were computed from a stale transform. A stopping branch returns its
    final classification state and has no next thinning pass.
    """
    if not isinstance(alignment_run, ToadsAtlasAlignmentRun):
        raise TypeError("alignment_run must be a selected registered TOADS run")
    if not isinstance(membership_pass, ToadsDefaultAtlasMembershipPass):
        raise TypeError("membership_pass must be a selected default-atlas membership pass")
    image_resolutions = tuple(
        float(value) for value in alignment_run.alignment.setup.signed_image_resolutions_mm
    )
    template_resolutions = tuple(
        float(value) for value in alignment_run.alignment.setup.template_resolutions_mm
    )
    atlas_shapes = toads_atlas_shapes_in_structure_order(
        alignment_run.atlas, alignment_run.atlas_data
    )
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    registration = continue_default_membership_through_registration_control(
        membership_pass,
        atlas_shapes,
        image_resolutions_mm=image_resolutions,
        template_resolutions_mm=template_resolutions,
        initial_transform=tuple(float(value) for value in transform),
        previous_energy=previous_energy,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        do_growing=do_growing,
        align_shapes=align_shapes,
        offset=offset,
        subsample=subsample,
        max_registration_iterations=max_registration_iterations,
        registration_energy_tolerance=registration_energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    registered = continue_default_membership_registered_iteration(
        registration,
        previous_energy=previous_energy,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volumes=class_volumes,
        boundary_distances=boundary_distances,
        template_labels=template_labels,
        previous_classification=previous_classification,
        relation_limit=relation_limit,
        do_growing=do_growing,
        initial_available=initial_available,
        record_growing_events=record_growing_events,
    )
    branch = registered.branch
    next_shapes = None
    next_thinning = None
    if not branch.control.stop:
        if branch.growing is None or branch.relations is None:
            raise RuntimeError("a continuing selected branch must grow and refresh relations")
        next_shapes = transform_toads_shape_priors(
            atlas_shapes,
            image_dimensions_xyz=tuple(int(value) for value in membership_pass.memberships.shape[:3]),
            image_resolutions_mm=image_resolutions,
            template_resolutions_mm=template_resolutions,
            transform=tuple(float(value) for value in registration.handoff.transform),
            class_bounds=registration.handoff.shape_bounds,
        )
        next_thinning = toads_post_registration_next_thinning(
            branch, next_shapes, record_thinning_items=record_thinning_items
        )

    stages = list(registration.stage_order)
    stages.extend(registered.stage_order[len(registration.stage_order):])
    if next_shapes is not None:
        stages.append(f"iteration{iteration + 1}.transformShapePriors")
    if next_thinning is not None:
        stages.extend((
            f"iteration{iteration + 1}.initThinningLabels",
            f"iteration{iteration + 1}.propagateThinning",
        ))
    final_segmentation = (
        branch.growing.segmentation.copy()
        if branch.growing is not None
        else np.asarray(segmentation).copy()
    )
    return ToadsRegisteredMembershipBranchCycle(
        membership_pass=membership_pass,
        registration_control=registration,
        registered_iteration=registered,
        next_transformed_shape_priors=next_shapes,
        next_thinning=next_thinning,
        final_segmentation=final_segmentation,
        stage_order=tuple(stages),
    )
