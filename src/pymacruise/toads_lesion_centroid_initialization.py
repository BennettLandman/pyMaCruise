"""Initialize the base-TOADS lesion centroid from selected atlas priors.

This ports the small wrapper stage that maps the per-channel lesion prior
through ``MultipleToadPreprocess.initialLesionCentroid``. The arithmetic is
shared with ordinary TOADS initial intensity centroids; this adapter makes the
distinct scalar-prior contract explicit.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_initialization import initial_intensity_centroids


def initial_toads_lesion_centroid(
    low: Any,
    high: Any,
    lesion_prior_by_channel: Any,
) -> NDArray[np.float32]:
    """Map one lesion prior per channel to the selected normalized intensity.

    The selected Java expression is ``Ilow[c] / (Ihigh[c] - Ilow[c]) +
    prior[c]``. ``low`` and ``high`` are the preprocessing robust bounds,
    computed before cropping, and ``lesion_prior_by_channel`` is the atlas
    lesion-prior vector in input-channel order. This routine does not use the
    image rotation or topology despite constructing a rotation in the Java
    method. As in ``initial_intensity_centroids``, non-finite inputs and
    nonpositive spans are rejected; zero-span Java NaN/Inf behavior is outside
    this bounded finite-input port.
    """
    low_values = np.asarray(low, dtype=np.float32)
    high_values = np.asarray(high, dtype=np.float32)
    prior_values = np.asarray(lesion_prior_by_channel, dtype=np.float32)
    if prior_values.ndim != 1:
        raise ValueError("lesion_prior_by_channel must be one-dimensional")
    if low_values.ndim != 1 or high_values.shape != low_values.shape:
        raise ValueError("low and high must be matching one-dimensional channel arrays")
    if prior_values.shape != low_values.shape:
        raise ValueError("lesion prior must have one value per image channel")
    return initial_intensity_centroids(
        low_values,
        high_values,
        prior_values[:, np.newaxis],
    )[:, 0].copy()
