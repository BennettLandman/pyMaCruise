"""Carry the source-faithful first TOADS update through the next thinning pass.

This adapter joins already validated statistics, membership, registration,
control, growth, relation, shape-resampling, and thinning operations. It is a
bounded registered-loop composition, not the complete TOADS workflow.
Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

from .toads_alignment_startup import ToadsAtlasAlignmentRun
from .toads_atlas import toads_atlas_shapes_in_structure_order
from .toads_default_membership_registration_control import (
    ToadsDefaultMembershipRegistrationControl,
    continue_default_membership_through_registration_control,
)
from .toads_default_membership_registered_iteration import (
    ToadsDefaultMembershipRegisteredIteration,
    continue_default_membership_registered_iteration,
)
from .toads_registered_first_exact_iteration import ToadsRegisteredFirstExactIteration
from .toads_shape_volume import transform_toads_shape_priors
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)


@dataclass(frozen=True)
class ToadsRegisteredFirstIterationCycle:
    """First exact memberships through registration, branch, and next thinning."""

    exact_iteration: ToadsRegisteredFirstExactIteration
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    next_thinning: ToadsPostRegistrationNextThinning
    next_transformed_shape_priors: NDArray[np.float32]
    stage_order: tuple[str, ...]


def continue_registered_first_exact_iteration_through_next_thinning(
    alignment_run: ToadsAtlasAlignmentRun,
    exact_iteration: ToadsRegisteredFirstExactIteration,
    *,
    maximum_iterations: int,
    convergence_threshold: float,
    relation_limit: int,
    do_growing: bool = False,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    record_growing_events: bool = True,
    record_thinning_items: bool = True,
    progress_callback: Callable[[str], None] | None = None,
) -> ToadsRegisteredFirstIterationCycle:
    """Run the first exact pass through the selected wrapper's next thinning.

    The first source loop starts at iteration one. Its previous energy is the
    initial membership objective, and its starting transform is the
    post-refresh atlas transform. Constructor-initialized grow ordering is
    zero. The branch uses the first thinning's label/skeleton/volume state and
    initialized relation distances; a continuing branch refreshes relations,
    transforms shape priors at the new rigid transform, and starts the next
    loop with thinning. Iteration and relation limits stay explicit because
    they are wrapper settings rather than image-derived values.
    """
    if not isinstance(alignment_run, ToadsAtlasAlignmentRun):
        raise TypeError("alignment_run must be a selected registered TOADS run")
    if not isinstance(exact_iteration, ToadsRegisteredFirstExactIteration):
        raise TypeError("exact_iteration must be the first exact registered TOADS iteration")
    if exact_iteration.startup is None:
        raise ValueError("exact_iteration must retain its registered startup state")

    startup = exact_iteration.startup
    membership = exact_iteration.membership_pass
    previous_energy = np.float32(
        startup.statistics_membership.membership_pass.objective
    )
    image_resolutions = tuple(
        float(value)
        for value in alignment_run.alignment.setup.signed_image_resolutions_mm
    )
    template_resolutions = tuple(
        float(value)
        for value in alignment_run.alignment.setup.template_resolutions_mm
    )
    atlas_shapes = toads_atlas_shapes_in_structure_order(
        alignment_run.atlas, alignment_run.atlas_data
    )
    initial_transform = tuple(
        float(value)
        for value in alignment_run.post_refresh_registration.final_transform
    )

    registration_control = continue_default_membership_through_registration_control(
        membership,
        atlas_shapes,
        image_resolutions_mm=image_resolutions,
        template_resolutions_mm=template_resolutions,
        initial_transform=initial_transform,
        previous_energy=previous_energy,
        iteration=1,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        do_growing=do_growing,
        align_shapes=align_shapes,
        offset=offset,
        subsample=subsample,
        max_registration_iterations=max_registration_iterations,
        registration_energy_tolerance=registration_energy_tolerance,
        minimum_lambda=minimum_lambda,
        progress_callback=progress_callback,
    )

    constructor = alignment_run.segmentation_constructor
    relation_state = alignment_run.relation_state
    initial_ordering = np.zeros(startup.thinning.segmentation.shape, dtype=np.float32)
    registered_iteration = continue_default_membership_registered_iteration(
        registration_control,
        previous_energy=previous_energy,
        iteration=1,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        segmentation=startup.thinning.segmentation,
        skeleton=startup.thinning.skeleton,
        ordering=initial_ordering,
        class_volumes=startup.thinning.class_volumes,
        boundary_distances=relation_state.boundary_distances,
        template_labels=np.asarray(
            [structure.label for structure in alignment_run.atlas.structures],
            dtype=np.int8,
        ),
        previous_classification=constructor.classification_labels,
        relation_limit=relation_limit,
        do_growing=do_growing,
        initial_available=startup.thinning.available,
        record_growing_events=record_growing_events,
        progress_callback=progress_callback,
    )
    branch = registered_iteration.branch
    if branch.control.stop:
        raise ValueError(
            "the first-iteration cycle ends at classification when the wrapper stops; "
            "use a later output adapter instead of requesting next thinning"
        )
    if branch.growing is None or branch.relations is None:
        raise RuntimeError("a continuing selected branch must grow and refresh relations")

    if progress_callback is not None:
        progress_callback("transforming TOADS shape priors after registration")
    next_shapes = transform_toads_shape_priors(
        atlas_shapes,
        image_dimensions_xyz=tuple(int(value) for value in membership.memberships.shape[:3]),
        image_resolutions_mm=image_resolutions,
        template_resolutions_mm=template_resolutions,
        transform=tuple(float(value) for value in registration_control.handoff.transform),
        class_bounds=registration_control.handoff.shape_bounds,
    )
    if progress_callback is not None:
        progress_callback("registered TOADS shape-prior transform complete")
        progress_callback("starting registered-cycle thinning")
    next_thinning = toads_post_registration_next_thinning(
        branch, next_shapes, record_thinning_items=record_thinning_items
    )
    if progress_callback is not None:
        progress_callback("registered-cycle thinning complete")
    stages = (
        *exact_iteration.stage_order,
        *registration_control.stage_order[1:],
        *registered_iteration.stage_order[len(registration_control.stage_order):],
        "iteration2.initThinningLabels",
        "iteration2.propagateThinning",
    )
    return ToadsRegisteredFirstIterationCycle(
        exact_iteration=exact_iteration,
        registration_control=registration_control,
        registered_iteration=registered_iteration,
        next_thinning=next_thinning,
        next_transformed_shape_priors=next_shapes,
        stage_order=stages,
    )
