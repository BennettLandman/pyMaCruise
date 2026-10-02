"""Run the selected registered default-atlas TOADS loop to its terminal branch.

This is orchestration over the individually checked selected stages; it does
not introduce replacement segmentation, registration, topology, or morphology
algorithms. The historical V3 BET/SPECTRE path is outside this native-RAS API.
Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_alignment_startup import ToadsAtlasAlignmentRun
from .toads_registered_first_exact_iteration import (
    ToadsRegisteredFirstExactIteration,
    continue_registered_unity_startup_through_first_exact_iteration,
)
from .toads_registered_first_iteration_cycle import (
    ToadsRegisteredFirstIterationCycle,
    continue_registered_first_exact_iteration_through_next_thinning,
)
from .toads_registered_membership_branch_cycle import (
    ToadsRegisteredMembershipBranchCycle,
    continue_registered_membership_through_branch,
)
from .toads_registered_following_exact_iteration import (
    continue_registered_cycle_through_following_exact_iteration,
)
from .toads_registered_unity_field_startup import (
    ToadsRegisteredUnityFieldStartup,
    compute_registered_toads_unity_field_startup,
)


@dataclass(frozen=True)
class ToadsRegisteredLoopConfig:
    """Explicit atlas/statistics parameters and selected loop controls."""

    category_exponents: Any
    exponent_by_channel_class: Any
    modality_weights: Any
    class_variance_weights: Any
    smoothing_priors: Any
    atlas_shape_priors: Any
    smoothing: float
    prior_factor: float
    prior_coefficient: float
    prior_scale: float = 0.05
    maximum_iterations: int = 50
    convergence_threshold: float = 1e-4
    relation_limit: int = 1
    do_growing: bool = False
    maximum_polynomial_degree: int = 3
    maximum_distance: float = 1e-4
    correct_inhomogeneity: bool = True
    centroid_smoothness: float = 0.1
    volume_ratio: float = 0.5
    thinning_limit: float = 3.0
    align_shapes: bool = True
    offset: int = 0
    subsample: int = 2
    max_registration_iterations: int = 5
    registration_energy_tolerance: float = 1e-6
    minimum_lambda: float = 1e-3


@dataclass(frozen=True)
class ToadsRegisteredLoopResult:
    """Source-ordered registered loop state through stop/classification."""

    alignment_run: ToadsAtlasAlignmentRun
    startup: ToadsRegisteredUnityFieldStartup
    exact_iterations: tuple[Any, ...]
    registered_cycles: tuple[Any, ...]
    classification: NDArray[np.int8] | None
    changed_classification_voxels: int | None


def run_registered_toads_loop(
    alignment_run: ToadsAtlasAlignmentRun,
    config: ToadsRegisteredLoopConfig,
) -> ToadsRegisteredLoopResult:
    """Execute registered TOADS startup and loop through the selected stop.

    The startup computes unity fields, initial exact memberships and first
    thinning. Each loop then carries the selected exact statistics/membership
    update through rigid registration and either next-loop thinning or final
    classification. ``maximum_iterations`` must be at least two because the
    current first-cycle adapter models its continuing branch explicitly.
    """
    if not isinstance(alignment_run, ToadsAtlasAlignmentRun):
        raise TypeError("alignment_run must be a selected registered TOADS run")
    if not isinstance(config, ToadsRegisteredLoopConfig):
        raise TypeError("config must be ToadsRegisteredLoopConfig")
    if (
        isinstance(config.maximum_iterations, bool)
        or not isinstance(config.maximum_iterations, int)
        or config.maximum_iterations < 2
    ):
        raise ValueError("maximum_iterations must be an integer of at least two")

    membership_options = {
        "exponent_by_channel_class": config.exponent_by_channel_class,
        "modality_weights": config.modality_weights,
        "smoothing_priors": config.smoothing_priors,
        "atlas_shape_priors": config.atlas_shape_priors,
        "smoothing": config.smoothing,
        "prior_factor": config.prior_factor,
        "prior_coefficient": config.prior_coefficient,
        "prior_scale": config.prior_scale,
    }
    startup = compute_registered_toads_unity_field_startup(
        alignment_run,
        **membership_options,
        volume_ratio=config.volume_ratio,
        thinning_limit=config.thinning_limit,
    )
    first_exact = continue_registered_unity_startup_through_first_exact_iteration(
        alignment_run,
        startup,
        category_exponents=config.category_exponents,
        class_variance_weights=config.class_variance_weights,
        centroid_smoothness=config.centroid_smoothness,
        **membership_options,
    )
    first_cycle = continue_registered_first_exact_iteration_through_next_thinning(
        alignment_run,
        first_exact,
        maximum_iterations=config.maximum_iterations,
        convergence_threshold=config.convergence_threshold,
        relation_limit=config.relation_limit,
        do_growing=config.do_growing,
        align_shapes=config.align_shapes,
        offset=config.offset,
        subsample=config.subsample,
        max_registration_iterations=config.max_registration_iterations,
        registration_energy_tolerance=config.registration_energy_tolerance,
        minimum_lambda=config.minimum_lambda,
    )
    exact_iterations: list[Any] = [first_exact]
    registered_cycles: list[Any] = [first_cycle]

    current_exact: Any = first_exact
    field_state = first_exact.fields
    current_cycle: ToadsRegisteredFirstIterationCycle | ToadsRegisteredMembershipBranchCycle = first_cycle
    current_iteration = 1
    while not current_cycle.registered_iteration.branch.control.stop:
        next_exact = continue_registered_cycle_through_following_exact_iteration(
            alignment_run,
            current_cycle,
            previous_degree=current_exact.degree,
            initial_centroids=current_exact.centroids,
            fields=field_state,
            category_exponents=config.category_exponents,
            class_variance_weights=config.class_variance_weights,
            maximum_polynomial_degree=config.maximum_polynomial_degree,
            maximum_distance=config.maximum_distance,
            correct_inhomogeneity=config.correct_inhomogeneity,
            centroid_smoothness=config.centroid_smoothness,
            **membership_options,
        )
        previous_energy = current_exact.membership_pass.objective
        next_iteration = current_iteration + 1
        previous_branch = current_cycle.registered_iteration.branch
        previous_thinning = current_cycle.next_thinning
        if previous_thinning is None or previous_branch.growing is None or previous_branch.relations is None:
            raise RuntimeError("a continuing registered branch must provide thinning and relation state")
        next_cycle = continue_registered_membership_through_branch(
            alignment_run,
            next_exact.membership_pass,
            previous_energy=previous_energy,
            iteration=next_iteration,
            maximum_iterations=config.maximum_iterations,
            convergence_threshold=config.convergence_threshold,
            segmentation=previous_thinning.thinning.segmentation,
            skeleton=previous_thinning.thinning.skeleton,
            ordering=previous_branch.growing.ordering,
            class_volumes=previous_thinning.thinning.class_volumes,
            boundary_distances=previous_branch.relations.boundary_distances,
            template_labels=np.asarray(
                [item.label for item in alignment_run.atlas.structures], dtype=np.int8
            ),
            previous_classification=alignment_run.segmentation_constructor.classification_labels,
            relation_limit=config.relation_limit,
            do_growing=config.do_growing,
            initial_transform=current_cycle.registration_control.handoff.transform,
            initial_available=previous_thinning.thinning.available,
            align_shapes=config.align_shapes,
            offset=config.offset,
            subsample=config.subsample,
            max_registration_iterations=config.max_registration_iterations,
            registration_energy_tolerance=config.registration_energy_tolerance,
            minimum_lambda=config.minimum_lambda,
        )
        exact_iterations.append(next_exact)
        registered_cycles.append(next_cycle)
        current_exact = next_exact
        current_cycle = next_cycle
        current_iteration = next_iteration

    final_branch = current_cycle.registered_iteration.branch
    return ToadsRegisteredLoopResult(
        alignment_run=alignment_run,
        startup=startup,
        exact_iterations=tuple(exact_iterations),
        registered_cycles=tuple(registered_cycles),
        classification=final_branch.classification,
        changed_classification_voxels=final_branch.changed_classification_voxels,
    )
