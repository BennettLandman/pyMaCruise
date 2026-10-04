"""Fit the selected exact polynomial correction field in channel order.

This bounded multichannel seam applies the existing selected single-channel
source kernel independently to each channel, carrying one field and centroid
vector per channel while sharing descriptor-ordered memberships and
segmentation. The 11-class base profile is the default; other class layouts
must provide their descriptor type names. It stops after field fitting;
statistics and membership accumulation remain separate validated stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_inhomogeneity_exact_fit import (
    ToadsExactChebyshevFieldFit,
    fit_toads_exact_chebyshev_field,
)


@dataclass(frozen=True)
class ToadsMultichannelExactFieldRefit:
    """Updated channel-last fields and selected fit outputs in input order."""

    fields: NDArray[np.float32]
    fits: tuple[ToadsExactChebyshevFieldFit, ...]
    channel_order: tuple[int, ...]


def toads_multichannel_exact_field_refit(
    image_values: Any,
    memberships: Any,
    segmentation: Any,
    centroids: Any,
    fields: Any,
    *,
    degree: int = 1,
    class_types: tuple[str, ...] | list[str] | None = None,
) -> ToadsMultichannelExactFieldRefit:
    """Refit one degree-1/2/3 field per channel in source channel order.

    The 11-class base atlas remains the default. Other atlas profiles must
    provide descriptor-ordered class types; this supports the selected
    13-class LesionTOADS correction profile without changing the base path.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    initial_fields = np.asarray(fields, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must be a channel-last 4-D array with at least one channel")
    channels = images.shape[-1]
    if mem.ndim != 4 or mem.shape[:3] != images.shape[:3] or mem.shape[-1] < 2:
        raise ValueError("memberships must match the image grid and contain at least two classes")
    classes = mem.shape[-1]
    if class_types is None and classes != 11:
        raise ValueError("class_types is required unless using the selected 11-class atlas")
    if class_types is not None and len(class_types) != classes:
        raise ValueError("class_types must contain one descriptor type per membership class")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be integer labels on the shared image grid")
    if centers.shape != (channels, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must have shape (channels, classes) in input channel order")
    if initial_fields.shape != images.shape or not np.isfinite(initial_fields).all():
        raise ValueError("fields must be finite and match channel-last image values")

    updated_fields = initial_fields.copy()
    fits: list[ToadsExactChebyshevFieldFit] = []
    order: list[int] = []
    for channel in range(channels):
        channel_field = updated_fields[..., channel].copy()
        fit = fit_toads_exact_chebyshev_field(
            images[..., channel],
            mem,
            labels,
            centers[channel],
            channel_field,
            degree=degree,
            class_types=class_types,
        )
        updated_fields[..., channel] = fit.field
        fits.append(fit)
        order.append(channel)
    return ToadsMultichannelExactFieldRefit(
        fields=updated_fields,
        fits=tuple(fits),
        channel_order=tuple(order),
    )
