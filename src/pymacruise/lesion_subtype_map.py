"""Selected LesionTOADS lesion subtype-map update.

This ports only ``LesionToadSegmentation.updateLesionClassifcation`` from the
selected VandyCRUISE source. The parent workflow invokes this map update after
ordinary labels are transferred, when its stopping condition is reached; this
module does not implement the enclosing segmentation workflow.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def update_lesion_subtype_map(
    classification: ArrayLike,
    template_labels: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    lesion_class: int,
    has_lesion: bool,
    lesion_classification: NDArray[np.int16],
) -> int:
    """Mutate the short subtype map and return the number of qualifying voxels.

    A voxel qualifies only when lesions are enabled, its current byte
    classification equals the configured lesion template label, and either
    lesion or black-hole membership is *strictly* greater than ordinary
    membership for that class. Qualifying voxels are subtype 2 only when
    black-hole membership is strictly greater than lesion membership; lesion
    wins exact subtype ties as subtype 1. Every nonqualifying voxel is cleared
    to zero, including stale values in the supplied output map.

    Volumes use the package's RAS-oriented ``(x, y, z)`` array order and the
    memberships array appends its class axis. Persisted image inputs and
    outputs should be loaded/saved through :mod:`pymacruise.image_io` as
    compressed ``.nii.gz`` images. The source returns a qualifying voxel count,
    not a count of map values changed.
    """
    labels = np.asarray(template_labels)
    classes = np.asarray(memberships, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)
    black_hole = np.asarray(black_hole_membership, dtype=np.float32)
    current = np.asarray(classification)

    if current.ndim != 3 or any(length == 0 for length in current.shape):
        raise ValueError("classification must be a nonempty 3-D (x, y, z) array")
    if current.dtype != np.dtype(np.int8):
        raise TypeError("classification must use signed int8 Java-byte values")
    if labels.ndim != 1 or not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("template_labels must be a one-dimensional integer array")
    if labels.size == 0 or np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template_labels must be nonempty signed Java-byte values")
    if classes.ndim != 4 or classes.shape[:3] != current.shape:
        raise ValueError("memberships must match classification with a class axis")
    if classes.shape[-1] != labels.size or classes.shape[-1] == 0:
        raise ValueError("memberships and template_labels must share a nonzero class count")
    if isinstance(lesion_class, (bool, np.bool_)) or not isinstance(
        lesion_class, (int, np.integer)
    ):
        raise TypeError("lesion_class must be an integer class index")
    if not 0 <= int(lesion_class) < classes.shape[-1]:
        raise ValueError("lesion_class is outside the membership class axis")
    if lesion.shape != current.shape or black_hole.shape != current.shape:
        raise ValueError("lesion and black-hole memberships must match classification")
    if not isinstance(has_lesion, (bool, np.bool_)):
        raise TypeError("has_lesion must be boolean")
    if not isinstance(lesion_classification, np.ndarray):
        raise TypeError("lesion_classification must be a writable numpy array")
    if lesion_classification.dtype != np.dtype(np.int16):
        raise TypeError("lesion_classification must use int16 Java-short storage")
    if lesion_classification.shape != current.shape:
        raise ValueError("lesion_classification must match classification")
    if not lesion_classification.flags.writeable:
        raise ValueError("lesion_classification must be writable")

    # Java's if/else loop clears the old map for all false predicates. Keep
    # comparisons separate: in particular, subtype ties favor lesion, while
    # qualification itself requires a strict win over ordinary membership.
    lesion_class_membership = classes[..., int(lesion_class)]
    qualifies = (
        bool(has_lesion)
        & (current == labels[int(lesion_class)])
        & (
            (lesion > lesion_class_membership)
            | (black_hole > lesion_class_membership)
        )
    )
    lesion_classification[...] = 0
    lesion_classification[qualifies] = 1
    black_hole_wins = qualifies & (lesion < black_hole)
    lesion_classification[black_hole_wins] = 2
    return int(np.count_nonzero(qualifies))
