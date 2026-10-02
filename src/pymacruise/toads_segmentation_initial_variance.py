"""Compute the selected TOADS startup class deviations after memberships."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_math import toads_optimized_class_standard_deviation


@dataclass(frozen=True)
class ToadsInitialClassVariances:
    """Per-channel, per-class standard deviations at wrapper startup."""

    standard_deviations: NDArray[np.float32]
    call_order: tuple[str, ...]


def initialize_toads_class_variances(
    *,
    images_xyz_channels: Any,
    memberships: Any,
    centroids: Any,
    segmentation: Any,
    class_types: Any,
) -> ToadsInitialClassVariances:
    """Run ``computeVariances`` for the selected non-lesion TOADS classes.

    The selected default atlas has MASK and optimized CSF/GM/WM classes. Their
    source exact-variance branches share the assigned-label formula. Startup
    correction fields, when enabled, are initialized to one, so this pass has
    the same arithmetic as the raw normalized image and avoids allocating a
    redundant full-grid field array.
    """
    images = np.asarray(images_xyz_channels, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    labels = np.asarray(segmentation)
    types = tuple(str(value).lower() for value in class_types)
    supported = {"mask", "optimized_csf", "optimized_gm", "optimized_wm"}
    if not types or any(value not in supported for value in types):
        raise ValueError("startup variance supports selected MASK and optimized tissue classes")
    if mems.ndim != 4 or mems.shape[-1] != len(types):
        raise ValueError("memberships must contain one class map per class type")
    deviations = toads_optimized_class_standard_deviation(
        images,
        mems,
        centers,
        labels,
        intensity_field=None,
    )
    return ToadsInitialClassVariances(
        deviations,
        ("segmentation.computeVariances", "segmentation.computeExactVariances"),
    )
