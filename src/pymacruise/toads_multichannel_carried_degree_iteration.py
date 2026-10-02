"""Carry TOADS exact-field degree through a multi-channel loop continuation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_default_exact_centroids import toads_default_exact_centroid_update
from .toads_default_exact_variances import toads_default_exact_class_variances
from .toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from .toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from .toads_multichannel_exact_field_refit import (
    ToadsMultichannelExactFieldRefit,
    toads_multichannel_exact_field_refit,
)
from .toads_post_registration_next_thinning import ToadsPostRegistrationNextThinning


@dataclass(frozen=True)
class ToadsMultichannelCarriedDegreeIteration:
    """One multi-channel exact field/statistics/membership continuation."""

    next_thinning: ToadsPostRegistrationNextThinning
    previous_degree: int
    degree: int
    energy_ratio: np.float32
    field_refit: ToadsMultichannelExactFieldRefit | None
    centroids: NDArray[np.float32]
    standard_deviations: NDArray[np.float32]
    image_variance_weights: NDArray[np.float32]
    image_standard_deviation: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def continue_toads_multichannel_through_carried_degree_iteration(
    next_thinning: ToadsPostRegistrationNextThinning,
    image_values: Any,
    initial_centroids: Any,
    fields: Any,
    *,
    previous_degree: int,
    polynomial_degree: int,
    low_energy_ratio: float,
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
) -> ToadsMultichannelCarriedDegreeIteration:
    """Continue a selected 11-class iteration with one field per image channel.

    Images and fields use ``(x,y,z,channel)``; centroids, exponents,
    modality weights, and class-variance weights use source channel order.
    A positive scheduled degree refits every channel, including an unchanged
    degree-three plateau. Degree zero leaves the attached fields unchanged.
    All atlas, relation, shape, and membership parameters remain explicit.
    """
    if not isinstance(next_thinning, ToadsPostRegistrationNextThinning):
        raise TypeError("next_thinning must be a selected post-registration thinning result")
    branch = next_thinning.branch
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("carried-degree iteration requires a continuing registered branch")

    images = np.asarray(image_values, dtype=np.float32)
    mem = np.asarray(branch.handoff.membership_pass.memberships, dtype=np.float32)
    labels = np.asarray(next_thinning.thinning.segmentation)
    centers0 = np.asarray(initial_centroids, dtype=np.float32)
    channel_fields = np.asarray(fields)
    categories = np.asarray(category_exponents, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    channels = images.shape[-1]
    if mem.shape != images.shape[:3] + (11,) or labels.shape != images.shape[:3]:
        raise ValueError("image, membership, and segmentation grids must match the selected 11-class layout")
    if centers0.shape != (channels, 11):
        raise ValueError("initial_centroids must have shape (channel,11)")
    if categories.shape != (channels, 3):
        raise ValueError("category_exponents must have shape (channel,CSF/GM/WM)")
    if channel_fields.shape != images.shape or channel_fields.dtype != np.float32:
        raise ValueError("fields must be a float32 array matching image_values")
    if not channel_fields.flags.writeable:
        raise ValueError("fields must be writable because the selected source updates them in place")

    ratio = np.float32(branch.control.energy_ratio)
    degree = toads_next_inhomogeneity_degree(
        previous_degree, polynomial_degree, ratio, low_energy_ratio
    )
    if degree not in (0, 1, 2, 3):
        raise ValueError("this selected default-atlas continuation supports degrees zero through three")

    field_refit = None
    stages = [*next_thinning.stage_order]
    if degree == 0:
        applied_fields = channel_fields
        stages.append("next_iteration.computeExact3DPolynomialField.degree0.noop")
    else:
        field_refit = toads_multichannel_exact_field_refit(
            images,
            mem,
            labels,
            centers0,
            channel_fields,
            degree=degree,
        )
        channel_fields[...] = field_refit.fields
        applied_fields = channel_fields
        stages.append(f"next_iteration.computeExact3DPolynomialField.degree{degree}")

    exponent_map = np.empty((channels, 11), dtype=np.float32)
    exponent_map[:, 0] = np.float32(1.0)
    exponent_map[:, 1:3] = categories[:, 0, None]
    exponent_map[:, 3:8] = categories[:, 1, None]
    exponent_map[:, 8:11] = categories[:, 2, None]
    supplied_exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    if supplied_exponents.shape != exponent_map.shape or not np.array_equal(supplied_exponents, exponent_map):
        raise ValueError("exponent_by_channel_class must match the default 11-class channel profile")

    centers = toads_default_exact_centroid_update(
        images,
        mem,
        centers0,
        categories,
        intensity_field=applied_fields,
        smoothness=centroid_smoothness,
    )
    deviations = toads_default_exact_class_variances(
        images,
        mem,
        centers,
        labels,
        intensity_field=applied_fields,
    )
    image_variance, image_std = toads_optimized_image_variance(
        images,
        mem[..., 1:],
        centers[:, 1:],
        intensity_field=applied_fields,
    )
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, image_variance, prior_scale
    )
    membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=mem,
        centroids=centers,
        exponent_by_channel_class=supplied_exponents,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variance,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=applied_fields,
    )
    stages.extend((
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    ))
    return ToadsMultichannelCarriedDegreeIteration(
        next_thinning=next_thinning,
        previous_degree=int(previous_degree),
        degree=degree,
        energy_ratio=ratio,
        field_refit=field_refit,
        centroids=centers,
        standard_deviations=deviations,
        image_variance_weights=image_variance,
        image_standard_deviation=image_std,
        centroid_similarity=similarity,
        membership_pass=membership,
        stage_order=tuple(stages),
    )
