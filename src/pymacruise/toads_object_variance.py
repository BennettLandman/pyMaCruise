"""Base TOADS assigned-class variance adapter for ordinary objects."""

from __future__ import annotations

from typing import TypeAlias

from numpy.typing import ArrayLike, NDArray
import numpy as np

from .toads_math import toads_optimized_class_standard_deviation


Array: TypeAlias = NDArray[np.float32]


def toads_object_assigned_class_standard_deviation(
    image_values: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    segmentation: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
) -> Array:
    """Compute exact assigned-class standard deviations for ordinary objects.

    This is a narrowly named adapter for the ordinary-object branch of
    ``MultipleToadSegmentation.computeExactVariances``. It delegates the
    arithmetic to the existing source-matched float32 implementation in
    :func:`toads_optimized_class_standard_deviation`; the optimized tissue
    assumptions of that helper are not needed here. Inputs must include only
    ordinary object classes, with arrays in ``(x, y, z, channel/class)`` order.

    The source scans the full grid, including boundary voxels assigned to the
    class, weights squared residuals by ``membership**2 + 1e-30f``, and returns
    zero for an unassigned class. When supplied, ``intensity_field`` is
    multiplied by the image before subtracting the centroid, in float32.
    Lesion object classes use additional lesion-membership routing and must not
    be passed to this adapter.
    """
    return toads_optimized_class_standard_deviation(
        image_values,
        memberships,
        centroids,
        segmentation,
        intensity_field=intensity_field,
    )
