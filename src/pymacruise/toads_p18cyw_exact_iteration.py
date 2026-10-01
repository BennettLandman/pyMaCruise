"""Continue P18cyv next-loop thinning through exact field/statistics/membership.

The selected AlgorithmMultipleToads order after propagateThinning is an
optional inhomogeneity-field update, exact centroids, exact variances, image
variance and centroid similarity during membership computation, then exact
memberships. This adapter composes those existing source-checked kernels with
fixed synthetic images/shape/closest-label inputs; it is not a complete TOADS
iteration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_p18cyv_next_iteration import ToadsP18cyvNextIteration
from .toads_post_thinning_exact_iteration import (
    ToadsPostThinningExactIteration,
    toads_post_thinning_exact_iteration_degree1,
)


@dataclass(frozen=True)
class ToadsP18cywExactIteration:
    """P18cyv thinning output continued through its selected exact iteration."""

    p18cyv: ToadsP18cyvNextIteration
    exact_iteration: ToadsPostThinningExactIteration
    stage_order: tuple[str, ...]


def continue_p18cyv_thinning_through_exact_iteration_degree1(
    p18cyv: ToadsP18cyvNextIteration,
    image_values: Any,
    initial_centroids: Any,
    field: Any,
    *,
    category_exponents: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    transformed_shapes: Any,
    closest_labels: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float,
    centroid_smoothness: float = 0.1,
) -> ToadsP18cywExactIteration:
    """Run the existing one-channel degree-one field/statistics/membership slice."""
    if not isinstance(p18cyv, ToadsP18cyvNextIteration):
        raise TypeError("p18cyv must be a selected P18cyv next-thinning result")
    next_thinning = p18cyv.next_thinning
    memberships = np.asarray(next_thinning.branch.handoff.membership_pass.memberships)
    image = np.asarray(image_values)
    if image.ndim != 3 or image.shape != memberships.shape[:3]:
        raise ValueError("image_values must be one channel on the P18cyv membership grid")
    labels = np.asarray(p18cyv.next_thinning.thinning.segmentation)
    if labels.shape != image.shape:
        raise ValueError("P18cyv thinning output must share the image grid")
    closest = np.asarray(closest_labels)
    if closest.shape != (*image.shape, 4):
        raise ValueError("closest_labels must contain the refreshed four-slot neighborhood map")

    exact = toads_post_thinning_exact_iteration_degree1(
        next_thinning,
        image_values,
        initial_centroids,
        field,
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=closest,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        centroid_smoothness=centroid_smoothness,
    )
    if exact.field_update.next_thinning is not next_thinning:
        raise RuntimeError("exact iteration did not preserve the P18cyv thinning input")
    return ToadsP18cywExactIteration(p18cyv, exact, exact.stage_order)
