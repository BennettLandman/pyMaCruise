"""Compute selected TOADS image-wide variance weights at initialization."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_math import toads_optimized_image_variance


@dataclass(frozen=True)
class ToadsInitialImageVariance:
    """Normalized inverse-variance channel weights and unnormalized stddev."""

    image_variance_weights: NDArray[np.float32]
    image_standard_deviations: NDArray[np.float32]
    call_order: tuple[str, ...]


def initialize_toads_image_variance(
    *,
    images_xyz_channels: Any,
    memberships: Any,
    centroids: Any,
    class_types: Any,
) -> ToadsInitialImageVariance:
    """Run the selected exact `computeImageVariance` tissue-class branch.

    The default descriptor has one MASK class followed by optimized CSF/GM/WM
    classes. The source image-variance loop excludes MASK and includes every
    optimized class over the full grid. Initial correction fields are all one,
    so this call uses the normalized source images directly.
    """
    images = np.asarray(images_xyz_channels, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    types = tuple(str(value).lower() for value in class_types)
    if len(types) < 2 or types[0] != "mask" or any(
        value not in {"optimized_csf", "optimized_gm", "optimized_wm"}
        for value in types[1:]
    ):
        raise ValueError("image variance expects the default MASK then optimized tissue class order")
    if mems.ndim != 4 or mems.shape[-1] != len(types):
        raise ValueError("memberships must contain one volume per class type")
    weights, standard_deviations = toads_optimized_image_variance(
        images,
        mems[..., 1:],
        centers[:, 1:],
        intensity_field=None,
    )
    return ToadsInitialImageVariance(
        weights,
        standard_deviations,
        (
            "segmentation.computeMemberships",
            "segmentation.computeImageVariance",
            "segmentation.computeExactImageVariance",
        ),
    )
