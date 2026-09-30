"""Image-wide variance update for selected base-TOADS ordinary OBJECT classes."""

from __future__ import annotations

from typing import Sequence

from numpy.typing import ArrayLike, NDArray
import numpy as np

from .toads_math import toads_optimized_image_variance


def toads_object_image_variance(
    image_values: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    *,
    class_types: Sequence[str],
    intensity_field: ArrayLike | None = None,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Return normalized ``varImg`` and image ``stddev`` for ordinary OBJECTs.

    This is the ordinary-OBJECT-only branch of selected VandyCRUISE
    ``MultipleToadSegmentation.computeExactImageVariance``. Arrays use
    ``(x, y, z, modality)`` images and ``(x, y, z, class)`` memberships.
    ``class_types`` is explicit so lesion and outlier mixture classes cannot
    silently enter this adapter. The source scans every voxel, including all
    image faces, accumulates weighted fourth-power residuals and denominators
    in float32, sets zero-numerator modality inverse variance to ``1e10f``,
    sums inverse variances in modality order, computes stddev from each raw
    inverse variance, and normalizes ``varImg`` by that sum.

    Optional ``intensity_field`` multiplies each image intensity before the
    class centroid residual, with float32 arithmetic. Computation is delegated
    to the already source-matched ``toads_optimized_image_variance`` kernel;
    the arithmetic is the same for selected optimized classes and ordinary
    OBJECTs. This helper does not implement lesion/outlier mixtures, approximate
    variance, or the surrounding membership estimation.

    Returns ``(normalized_inverse_variance, image_stddev)``; stddev is
    calculated before normalization, matching the Java method.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    if images.ndim != 4 or mems.ndim != 4:
        raise ValueError("images and memberships must be 4-D arrays")
    if mems.shape[:3] != images.shape[:3]:
        raise ValueError("memberships and images must share a spatial grid")
    if len(class_types) != mems.shape[-1]:
        raise ValueError("class_types must name every membership class")
    if any(str(kind).upper() != "OBJECT" for kind in class_types):
        raise ValueError("only ordinary OBJECT classes are supported")
    return toads_optimized_image_variance(
        images,
        mems,
        centers,
        intensity_field=intensity_field,
    )
