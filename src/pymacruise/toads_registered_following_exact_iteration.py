"""Continue a registered loop cycle through its following exact update."""

from __future__ import annotations

from typing import Any

import numpy as np

from .toads_alignment_startup import ToadsAtlasAlignmentRun
from .toads_registered_first_iteration_cycle import ToadsRegisteredFirstIterationCycle
from .toads_multichannel_carried_degree_iteration import (
    ToadsMultichannelCarriedDegreeIteration,
    continue_toads_multichannel_through_carried_degree_iteration,
)
from .toads_registered_membership_branch_cycle import ToadsRegisteredMembershipBranchCycle


def continue_registered_cycle_through_following_exact_iteration(
    alignment_run: ToadsAtlasAlignmentRun,
    registered_cycle: ToadsRegisteredMembershipBranchCycle,
    *,
    previous_degree: int,
    initial_centroids: Any,
    fields: Any,
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
    """Apply the next thinning's degree decision and exact statistics pass.

    A stopped registered cycle has no following loop pass and is rejected.
    On a continuing branch, the existing carried-degree adapter performs the
    selected field no-op/refit followed by centroid, variance, similarity, and
    membership updates. The energy ratio and transformed spatial priors are
    those from the immediately preceding registered cycle.
    """
    if not isinstance(alignment_run, ToadsAtlasAlignmentRun):
        raise TypeError("alignment_run must be a selected registered TOADS run")
    if not isinstance(
        registered_cycle,
        (ToadsRegisteredFirstIterationCycle, ToadsRegisteredMembershipBranchCycle),
    ):
        raise TypeError("registered_cycle must be a selected continuing registered cycle")
    branch = registered_cycle.registered_iteration.branch
    if (
        registered_cycle.next_thinning is None
        or branch.control.stop
        or branch.relations is None
    ):
        raise ValueError("the selected registered branch stopped before another exact iteration")
    images = np.moveaxis(
        np.asarray(alignment_run.template_refresh.recropped_images, dtype=np.float32),
        0,
        -1,
    )
    if images.shape != np.asarray(fields).shape:
        raise ValueError("registered images and attached fields must share channel-last dimensions")
    maximum_distance32 = np.float32(maximum_distance)
    if isinstance(maximum_distance, (bool, np.bool_)) or not np.isfinite(maximum_distance32) or maximum_distance32 < 0:
        raise ValueError("maximum_distance must be finite and nonnegative")
    low_energy_ratio = float(np.sqrt(maximum_distance32))
    return continue_toads_multichannel_through_carried_degree_iteration(
        registered_cycle.next_thinning,
        images,
        initial_centroids,
        fields,
        previous_degree=previous_degree,
        polynomial_degree=maximum_polynomial_degree,
        low_energy_ratio=low_energy_ratio,
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        transformed_shapes=registered_cycle.next_transformed_shape_priors,
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
