"""Continue the registered first-loop cycle through iteration-two memberships."""

from __future__ import annotations

from typing import Any

import numpy as np

from .toads_alignment_startup import ToadsAtlasAlignmentRun
from .toads_multichannel_carried_degree_iteration import (
    ToadsMultichannelCarriedDegreeIteration,
    continue_toads_multichannel_through_carried_degree_iteration,
)
from .toads_registered_first_iteration_cycle import ToadsRegisteredFirstIterationCycle


def continue_registered_first_cycle_through_second_exact_iteration(
    alignment_run: ToadsAtlasAlignmentRun,
    first_cycle: ToadsRegisteredFirstIterationCycle,
    *,
    category_exponents: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float,
    maximum_polynomial_degree: int = 3,
    maximum_distance: float = 1e-4,
    correct_inhomogeneity: bool = True,
    centroid_smoothness: float = 0.1,
) -> ToadsMultichannelCarriedDegreeIteration:
    """Run iteration-two field scheduling, statistics, and exact memberships.

    This maps the P18cybf next-thinning state into the existing selected
    multichannel continuation. The wrapper's inhomogeneity threshold is
    ``sqrt(maximum_distance)``; degree is carried from the first exact pass
    and updated only when correction is enabled. All statistical and prior
    settings remain explicit because they are selected wrapper parameters.
    """
    if not isinstance(alignment_run, ToadsAtlasAlignmentRun):
        raise TypeError("alignment_run must be a selected registered TOADS run")
    if not isinstance(first_cycle, ToadsRegisteredFirstIterationCycle):
        raise TypeError("first_cycle must be a completed registered first-iteration cycle")
    branch = first_cycle.registered_iteration.branch
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("iteration two requires the first cycle's continuing registered branch")
    images = np.moveaxis(
        np.asarray(alignment_run.template_refresh.recropped_images, dtype=np.float32),
        0,
        -1,
    )
    if images.shape != first_cycle.exact_iteration.fields.shape:
        raise ValueError("registered images and carried correction fields must share channel-last dimensions")
    if isinstance(maximum_distance, (bool, np.bool_)):
        raise ValueError("maximum_distance must be finite and nonnegative")
    maximum_distance32 = np.float32(maximum_distance)
    if not np.isfinite(maximum_distance32) or maximum_distance32 < 0:
        raise ValueError("maximum_distance must be finite and nonnegative")
    low_energy_ratio = float(np.sqrt(maximum_distance32))
    return continue_toads_multichannel_through_carried_degree_iteration(
        first_cycle.next_thinning,
        images,
        first_cycle.exact_iteration.centroids,
        first_cycle.exact_iteration.fields,
        previous_degree=first_cycle.exact_iteration.degree,
        polynomial_degree=maximum_polynomial_degree,
        low_energy_ratio=low_energy_ratio,
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        transformed_shapes=first_cycle.next_transformed_shape_priors,
        closest_labels=branch.relations.closest_labels,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        centroid_smoothness=centroid_smoothness,
        correct_inhomogeneity=correct_inhomogeneity,
    )
