"""Source-faithful ordinary LesionTOADS membership classification update.

The selected VandyCRUISE ``LesionToadSegmentation.updateClassificationFromMemberships``
method mutates its byte classification volume in place. This adapter covers only
that method's voxelwise class selection; it is not the LesionTOADS algorithm.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def update_lesion_classification_from_memberships(
    classification: NDArray[np.int8],
    template_labels: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    lesion_class: int,
) -> int:
    """Update an ordinary byte classification volume in place and return changes.

    At every ``(x, y, z)`` voxel Java starts with ordinary class 0 and scans
    classes 1 onward using strict ``>`` comparisons, so the first class wins
    exact ties. Lesion wins only when its membership is strictly greater than
    the best ordinary membership. Otherwise, any positive winning ordinary
    membership is mapped through ``template_labels``, including class index 0;
    a nonpositive winner produces label zero. The return value is the number of byte values changed
    from the input classification, including changes to/from signed labels.

    ``classification``, ``memberships`` and ``lesion_membership`` use the
    package's RAS-oriented ``(x, y, z)`` array convention; memberships append a
    class axis. ``classification`` must be a writable ``numpy.int8`` array so
    that mutation and Java signed-byte behavior remain explicit. The arrays
    contain voxel values only; load/save NIfTI through ``image_io`` to enforce
    RAS orientation and compressed ``.nii.gz`` storage.
    """
    if not isinstance(classification, np.ndarray):
        raise TypeError("classification must be a writable numpy.int8 array")
    if classification.dtype != np.dtype(np.int8):
        raise TypeError("classification must use signed int8 Java-byte storage")
    if not classification.flags.writeable:
        raise ValueError("classification must be writable for the in-place update")
    if classification.ndim != 3 or any(size == 0 for size in classification.shape):
        raise ValueError("classification must be a nonempty 3-D (x, y, z) array")

    labels = np.asarray(template_labels)
    if labels.ndim != 1 or not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("template_labels must be a one-dimensional integer array")
    if labels.size == 0:
        raise ValueError("template_labels must contain at least one class")
    if np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template_labels must fit Java's signed byte range")

    ordinary = np.asarray(memberships, dtype=np.float32)
    if ordinary.ndim != 4 or ordinary.shape[:3] != classification.shape:
        raise ValueError("memberships must match classification with a class axis")
    classes = ordinary.shape[-1]
    if classes == 0 or classes != labels.size:
        raise ValueError("memberships and template_labels must have the same nonzero class count")
    if isinstance(lesion_class, (bool, np.bool_)) or not isinstance(
        lesion_class, (int, np.integer)
    ):
        raise TypeError("lesion_class must be an integer class index")
    if not 0 <= int(lesion_class) < classes:
        raise ValueError("lesion_class is outside the membership class axis")

    lesion = np.asarray(lesion_membership, dtype=np.float32)
    if lesion.shape != classification.shape:
        raise ValueError("lesion_membership must match classification")
    if not np.isfinite(ordinary).all() or not np.isfinite(lesion).all():
        raise ValueError("memberships must contain only finite values")

    previous = classification.copy()
    # numpy.argmax returns the first maximum, matching Java's strict-`>` scan.
    best = np.argmax(ordinary, axis=-1)
    best_membership = np.take_along_axis(ordinary, best[..., None], axis=-1)[..., 0]
    updated = np.zeros(classification.shape, dtype=np.int8)

    lesion_wins = lesion > best_membership
    updated[lesion_wins] = np.int8(labels[int(lesion_class)])

    ordinary_wins = (~lesion_wins) & (best_membership > 0)
    if np.any(ordinary_wins):
        updated[ordinary_wins] = labels[best[ordinary_wins]].astype(np.int8, copy=False)

    classification[...] = updated
    return int(np.count_nonzero(previous != classification))
