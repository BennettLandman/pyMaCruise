"""Continue P18czj thinning through the next carried degree-one iteration.

The selected field/statistics/membership kernels are reused. The prior
correction field, exact centroids, and polynomial degree are carried from
P18czh across registration, branch refresh, and thinning.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .toads_p18cza_carried_degree_iteration import (
    ToadsP18czaCarriedDegreeIteration,
    continue_p18cyz_thinning_through_carried_degree1_iteration,
)
from .toads_p18czj_next_loop_thinning import ToadsP18czjNextLoopThinning


@dataclass(frozen=True)
class ToadsP18czkNextExactIteration:
    """P18czj next thinning attached to its following exact iteration."""

    previous_loop: ToadsP18czjNextLoopThinning
    exact_iteration: ToadsP18czaCarriedDegreeIteration
    stage_order: tuple[str, ...]


def continue_p18czj_thinning_through_next_exact_iteration_degree1(
    previous_loop: ToadsP18czjNextLoopThinning,
    image_values: Any,
    *,
    polynomial_degree: int,
    low_energy_ratio: float,
    category_exponents: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    transformed_shapes: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float,
    centroid_smoothness: float = 0.1,
) -> ToadsP18czkNextExactIteration:
    """Run source-ordered degree-one field, statistics, and memberships."""
    if not isinstance(previous_loop, ToadsP18czjNextLoopThinning):
        raise TypeError("previous_loop must be a selected P18czj result")

    registered = previous_loop.registered_branch
    previous_exact = registered.exact_iteration.exact_iteration
    previous_stats = previous_exact.exact_iteration
    branch = previous_loop.next_thinning.branch
    if branch.control.stop:
        raise ValueError("next exact iteration requires a continuing P18czi branch")
    if branch.growing is None or branch.relations is None:
        raise ValueError("next exact iteration requires carried growth and refreshed relations")

    exact = continue_p18cyz_thinning_through_carried_degree1_iteration(
        previous_loop.next_thinning,
        image_values,
        previous_stats.centroids[0],
        previous_stats.field_update.fit.field,
        previous_degree=previous_exact.degree,
        polynomial_degree=polynomial_degree,
        low_energy_ratio=low_energy_ratio,
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=branch.relations.closest_labels,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        centroid_smoothness=centroid_smoothness,
    )
    return ToadsP18czkNextExactIteration(previous_loop, exact, exact.stage_order)
