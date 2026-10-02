"""Carry selected atlas-startup outputs into the first exact-membership pass.

This adapter connects the registered, recropped RAS state produced by
``register_toads_from_atlas_assets`` to the existing source-ordered startup
statistics and membership composition. Atlas-derived arrays are taken from
that result; numerical atlas/segmentation settings remain explicit so this
does not invent defaults or imply that field fitting and TOADS iteration have
been assembled.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_initial_statistics_membership import (
    ToadsInitialStatisticsMembership,
    compute_toads_initial_statistics_membership,
)


def compute_registered_toads_initial_membership(
    alignment_run: Any,
    *,
    exponent_by_channel_class: ArrayLike,
    modality_weights: ArrayLike,
    smoothing_priors: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float = 0.05,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsInitialStatisticsMembership:
    """Run initial statistics and exact membership from registered atlas state.

    ``alignment_run`` must have been produced by
    :func:`register_toads_from_atlas_assets` using the selected default
    11-class descriptor. Its channel-first recropped image is converted to
    RAS ``(x,y,z,channel)`` order. Memberships, centroids, constructor labels,
    descriptor class names, relation-derived closest labels, and transformed
    shape-prior samples are taken from the run object. Exponents, tissue
    modality weights, atlas shape/smoothing weights, and effective
    segmentation coefficients stay explicit inputs because they are
    descriptor or wrapper settings, not quantities inferred by registration.

    This closes the atlas-startup-to-membership data handoff only. It does
    not estimate correction fields, run an outer TOADS iteration, or claim a
    combined Java wrapper comparison.
    """
    try:
        atlas = alignment_run.atlas
        structures = atlas.structures
        recropped = np.asarray(
            alignment_run.template_refresh.recropped_images, dtype=np.float32
        )
        memberships = np.asarray(
            alignment_run.segmentation_initialization.memberships, dtype=np.float32
        )
        centroids = np.asarray(
            alignment_run.segmentation_initialization.centroids, dtype=np.float32
        )
        segmentation = np.asarray(
            alignment_run.segmentation_constructor.segmentation
        )
        closest_labels = np.asarray(alignment_run.relation_state.closest_labels)
        transformed_shapes = np.asarray(
            alignment_run.transformed_shape_priors, dtype=np.float32
        )
    except AttributeError as exc:
        raise TypeError(
            "alignment_run must contain atlas, template_refresh, "
            "segmentation_initialization, segmentation_constructor, "
            "relation_state, and transformed_shape_priors state"
        ) from exc

    if recropped.ndim != 4 or recropped.shape[0] < 1:
        raise ValueError("alignment_run recropped images must have shape (channels,x,y,z)")
    class_types = tuple(str(structure.topology).upper() for structure in structures)
    if class_types != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError(
            "registered initial exact membership supports the selected default 11-class atlas order"
        )
    image_values = np.moveaxis(recropped, 0, -1)
    if memberships.shape != image_values.shape[:3] + (len(class_types),):
        raise ValueError("registered memberships must match the recropped image and atlas grid")
    if centroids.shape != (image_values.shape[-1], len(class_types)):
        raise ValueError("registered centroids must match the image channels and atlas classes")
    if segmentation.shape != image_values.shape[:3]:
        raise ValueError("constructor segmentation must match the recropped image grid")
    if transformed_shapes.shape != memberships.shape:
        raise ValueError("transformed shape priors must match registered memberships")
    if closest_labels.shape != image_values.shape[:3] + (4,):
        raise ValueError("relation closest-label map must match the recropped image grid")

    return compute_toads_initial_statistics_membership(
        image_values=image_values,
        initial_memberships=memberships,
        centroids=centroids,
        segmentation=segmentation,
        class_types=class_types,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
