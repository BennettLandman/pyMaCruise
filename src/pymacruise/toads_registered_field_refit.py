"""Fit selected TOADS correction fields from registered startup state."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_inhomogeneity_exact_fit import ToadsExactChebyshevFieldFit
from .toads_multichannel_exact_field_refit import toads_multichannel_exact_field_refit


@dataclass(frozen=True)
class ToadsRegisteredFieldRefit:
    """Source-ordered field fits using an atlas-registration startup result."""

    fields: NDArray[np.float32]
    fits: tuple[ToadsExactChebyshevFieldFit, ...]
    degree: int
    channel_order: tuple[int, ...]
    stage_order: tuple[str, ...]


def fit_registered_toads_correction_fields(
    alignment_run: Any,
    *,
    degree: int = 3,
) -> ToadsRegisteredFieldRefit:
    """Fit the exact global correction field from selected startup state.

    The selected VandyCRUISE base wrapper constructs a degree-three global
    Chebyshev correction after post-registration image/segmentation
    initialization. This adapter consumes the corresponding recropped RAS
    images, exact membership maps, constructor labels, and centroids from one
    ``ToadsAtlasAlignmentRun`` and starts every attached field at unity.
    Degree and every numerical setting remain in the source-compatible field
    kernel; this function only carries the wrapper state and channel order.
    """
    if isinstance(degree, bool) or not isinstance(degree, (int, np.integer)) or int(degree) not in (1, 2, 3):
        raise ValueError("degree must be one of the selected supported values: 1, 2, or 3")
    degree = int(degree)
    try:
        class_types = tuple(str(item.topology).upper() for item in alignment_run.atlas.structures)
        images_channel_first = np.asarray(
            alignment_run.template_refresh.recropped_images, dtype=np.float32
        )
        memberships = np.asarray(
            alignment_run.segmentation_initialization.memberships, dtype=np.float32
        )
        segmentation = np.asarray(alignment_run.segmentation_constructor.segmentation)
        centroids = np.asarray(
            alignment_run.segmentation_initialization.centroids, dtype=np.float32
        )
    except AttributeError as exc:
        raise TypeError("alignment_run is missing selected registered TOADS startup state") from exc

    if class_types != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError("registered correction fields support the selected default 11-class atlas order")
    if images_channel_first.ndim != 4 or images_channel_first.shape[0] < 1:
        raise ValueError("registered recropped images must have shape (channels,x,y,z)")
    images = np.moveaxis(images_channel_first, 0, -1)
    if memberships.shape != images.shape[:3] + (11,):
        raise ValueError("registered memberships must match the recropped image and 11-class atlas")
    if segmentation.shape != images.shape[:3] or not np.issubdtype(segmentation.dtype, np.integer):
        raise ValueError("constructor segmentation must be integer labels on the registered image grid")
    if centroids.shape != (images.shape[-1], 11):
        raise ValueError("registered centroids must match channel count and atlas class order")

    initial_fields = np.ones(images.shape, dtype=np.float32)
    refit = toads_multichannel_exact_field_refit(
        images,
        memberships,
        segmentation,
        centroids,
        initial_fields,
        degree=degree,
    )
    return ToadsRegisteredFieldRefit(
        fields=refit.fields,
        fits=refit.fits,
        degree=degree,
        channel_order=refit.channel_order,
        stage_order=(f"registered_startup.computeExact3DPolynomialField.degree{degree}",),
    )
