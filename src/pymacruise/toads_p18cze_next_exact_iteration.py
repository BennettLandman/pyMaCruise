"""Continue P18czd thinning through the next carried degree-one iteration.

This orchestrates existing exact-field, statistics, and membership kernels.
The correction degree and prior field/centroids are taken from the immediately
preceding P18cza iteration; class priors and transformed shapes remain explicit
fixed synthetic inputs. This is one selected single-channel loop segment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .toads_p18cza_carried_degree_iteration import (
    ToadsP18czaCarriedDegreeIteration,
    continue_p18cyz_thinning_through_carried_degree1_iteration,
)
from .toads_p18czd_next_loop_thinning import ToadsP18czdNextLoopThinning


@dataclass(frozen=True)
class ToadsP18czeNextExactIteration:
    """Next exact iteration attached to its P18czd next-loop thinning input."""

    previous_thinning: ToadsP18czdNextLoopThinning
    exact_iteration: ToadsP18czaCarriedDegreeIteration
    stage_order: tuple[str, ...]


def continue_p18czd_thinning_through_next_exact_iteration_degree1(
    previous_thinning: ToadsP18czdNextLoopThinning,
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
) -> ToadsP18czeNextExactIteration:
    """Run the carried degree gate, field refit, statistics, and memberships.

    The incoming degree, attached field, and initial centroids are read from
    the P18cza exact output that has been carried through registration,
    branch refresh, and P18czd thinning. The energy ratio comes from that
    intervening branch's membership control state, as in the source loop.
    """
    if not isinstance(previous_thinning, ToadsP18czdNextLoopThinning):
        raise TypeError("previous_thinning must be a selected P18czd result")
    registered = previous_thinning.registered_branch
    previous_exact = registered.exact_iteration
    previous_stats = previous_exact.exact_iteration
    branch = previous_thinning.next_thinning.branch
    if branch.control.stop:
        raise ValueError("next-loop exact iteration requires a continuing branch")
    if branch.relations is None or branch.growing is None:
        raise ValueError("next-loop exact iteration requires refreshed relation/growth state")

    exact = continue_p18cyz_thinning_through_carried_degree1_iteration(
        previous_thinning.next_thinning,
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
    return ToadsP18czeNextExactIteration(
        previous_thinning=previous_thinning,
        exact_iteration=exact,
        stage_order=exact.stage_order,
    )
