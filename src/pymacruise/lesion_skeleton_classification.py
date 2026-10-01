"""Selected LesionTOADS skeleton-to-label classification update.

Shiee et al. (2010), ``docs/papers/lesion-toads.md``, discuss preserving
topology in lesion-aware segmentation. This skeleton-to-label byte update is
only a selected VandyCRUISE output helper; its sentinel and class-index
semantics are source-specific and are not specified by the paper. The helper
is dormant in the selected wrapper, and its synthetic Java-excerpt validation
does not show paper-level topology preservation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def update_lesion_classification_from_skeleton(
    classification: NDArray[np.int8],
    template_labels: ArrayLike,
    skeleton: ArrayLike,
) -> int:
    """Overwrite classification from a Java-byte skeleton and return changes.

    The selected source uses ``EMPTY=-1``. Empty skeleton voxels become byte
    zero; all other skeleton values index ``template_labels``. The return is
    the number of voxels whose signed-byte label changed. Arrays use the
    package's ``(x, y, z)`` convention; RAS NIfTI I/O stays in ``image_io``.

    This ports the method only. It does not claim that the selected workflow
    calls this dormant helper.
    """
    if not isinstance(classification, np.ndarray) or classification.dtype != np.int8:
        raise TypeError("classification must be a writable numpy.int8 array")
    if not classification.flags.writeable:
        raise ValueError("classification must be writable")
    if classification.ndim != 3 or any(n == 0 for n in classification.shape):
        raise ValueError("classification must be a nonempty 3-D array")

    labels = np.asarray(template_labels)
    if labels.ndim != 1 or labels.size == 0 or not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("template_labels must be a nonempty one-dimensional integer array")
    if np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template_labels must fit Java's signed byte range")

    skel = np.asarray(skeleton)
    if skel.dtype != np.int8:
        raise TypeError("skeleton must use signed int8 Java-byte storage")
    if skel.shape != classification.shape:
        raise ValueError("skeleton must match classification shape")
    nonempty = skel != np.int8(-1)
    if np.any(skel[nonempty] < 0) or np.any(skel[nonempty] >= labels.size):
        raise ValueError("nonempty skeleton values must index template_labels")

    before = classification.copy()
    classification[...] = 0
    if np.any(nonempty):
        classification[nonempty] = labels.astype(np.int8, copy=False)[skel[nonempty]]
    return int(np.count_nonzero(before != classification))
