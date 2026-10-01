"""A bounded, source-ordered continuing TOADS iteration slice.

This adapter composes already validated kernels for thinning, one ordinary
OBJECT membership sweep, the non-stop decision, populated growing, and relation
refresh. It is an integrated slice only: the upstream centroid/variance
dispatcher, atlas registration, optional inhomogeneity correction, and full
TOADS initialization/convergence workflow are deliberately excluded. Model
parameters are supplied by the caller, as if produced by the preceding
source stages.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
an iterative combination of statistical classification and topology-aware
updates. This adapter composes one continuing source-defined pass in that
general role. The tested control, growth, and relation-refresh sequence is
selected VandyCRUISE behavior with fixed model parameters; it is not a full
TOADS iteration from image initialization through convergence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_control import ToadsIterationControl, toads_iteration_control
from pymacruise.toads_math import toads_object_membership_pass
from pymacruise.toads_relations import (
    ToadsRelationState,
    toads_refresh_relation_state,
)
from pymacruise.toads_thinning import (
    ToadsGrowingPass,
    ToadsThinningPass,
    toads_propagate_growing,
    toads_propagate_thinning,
)


@dataclass(frozen=True)
class ToadsContinuingIterationSlice:
    """Fieldwise outputs of the bounded continuing-iteration composition."""

    thinning: ToadsThinningPass
    memberships: NDArray[np.float32]
    membership_energy: np.float32
    control: ToadsIterationControl
    growing: ToadsGrowingPass
    relations: ToadsRelationState


def toads_continuing_iteration_slice(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    relation_matrix: ArrayLike,
    thinning_shape_values: ArrayLike,
    volume_ratio: float,
    thinning_limit: float,
    image_values: ArrayLike,
    centroids: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    membership_shape_values: ArrayLike,
    membership_closest_labels: ArrayLike,
    smoothing_priors: ArrayLike,
    centroid_similarity: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    ordering: ArrayLike,
    boundary_distances: ArrayLike,
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    *,
    relation_limit: int,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsContinuingIterationSlice:
    """Run a supported non-stopping iteration over a supplied 3-D state.

    Source order is ``propagateThinning``; ordinary OBJECT membership sweep;
    source energy/iteration decision; ``propagateGrowing``; then relations,
    boundary distances, and closest labels. The caller must provide the
    relation/shape model for thinning, the ordinary OBJECT membership model,
    and the preceding iteration's ordering and boundary-distance fields.
    This function refuses a stopping decision because classification transfer
    belongs to the separate stopping branch.
    """
    thin = toads_propagate_thinning(
        segmentation,
        memberships,
        relation_matrix,
        thinning_shape_values,
        volume_ratio,
        thinning_limit,
    )
    updated_memberships, objective = toads_object_membership_pass(
        image_values,
        memberships,
        centroids,
        modality_weights,
        class_variance_weights,
        image_variances,
        membership_shape_values,
        membership_closest_labels,
        smoothing_priors,
        centroid_similarity,
        atlas_shape_priors,
        smoothing,
        prior_factor,
        prior_coefficient,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    # The source compares the previous energy with the updated sweep objective.
    # Growing runs on every non-stop iteration, even when doGrowing is false.
    control = toads_iteration_control(
        previous_energy,
        objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=False,
    )
    if control.stop:
        raise ValueError("continuing iteration slice requires a non-stopping decision")
    growing = toads_propagate_growing(
        thin.segmentation,
        thin.skeleton,
        updated_memberships,
        ordering,
        thin.class_volumes,
        initial_available=thin.available,
    )
    relations = toads_refresh_relation_state(
        growing.segmentation,
        updated_memberships,
        boundary_distances,
        relation_limit,
    )
    return ToadsContinuingIterationSlice(
        thinning=thin,
        memberships=updated_memberships,
        membership_energy=objective,
        control=control,
        growing=growing,
        relations=relations,
    )
