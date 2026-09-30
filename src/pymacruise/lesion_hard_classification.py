"""Selected LesionTOADS membership-to-label output map.

This ports the voxel selection in ``LesionToadSegmentation.generateClassification``.
It maps ordinary memberships to atlas labels, then lets the lesion-membership
map override only when it is strictly greater than the best ordinary score.
Workflow classification updates and NIfTI serialization are separate stages.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def generate_lesion_toads_classification(
    template_labels: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    lesion_class: int,
) -> NDArray[np.int8]:
    """Return the selected Java method's signed-byte classification map.

    ``memberships`` has shape ``(x, y, z, class)`` and ``lesion_membership``
    has shape ``(x, y, z)``. The ordinary scan starts at score zero and updates
    only for strict ``>`` comparisons, so the first positive class wins ties,
    and zero/negative/NaN-only scores leave no ordinary winner. Lesion wins
    only when its float32 score is strictly greater than the best ordinary
    score. Equality therefore preserves the ordinary label. Positive
    infinity and NaN follow Java float comparison behavior.

    ``template_labels`` stores signed Java byte labels, with ``lesion_class``
    identifying the label used when lesion membership wins. Arrays are voxel
    values in the package's ``(x, y, z)`` convention; NIfTI I/O and RAS
    orientation are handled by :mod:`pymacruise.image_io`.
    """
    labels = np.asarray(template_labels)
    ordinary = np.asarray(memberships, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)

    if labels.ndim != 1 or not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("template_labels must be a one-dimensional integer array")
    if labels.size == 0 or np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template_labels must contain signed Java byte labels")
    if ordinary.ndim != 4 or any(n == 0 for n in ordinary.shape):
        raise ValueError("memberships must be a nonempty (x, y, z, class) array")
    if ordinary.shape[-1] == 0 or ordinary.shape[-1] > labels.size:
        raise ValueError("membership classes must index template_labels")
    if lesion.shape != ordinary.shape[:3]:
        raise ValueError("lesion_membership must match the membership spatial shape")
    if isinstance(lesion_class, (bool, np.bool_)) or not isinstance(
        lesion_class, (int, np.integer)
    ):
        raise TypeError("lesion_class must be an integer label index")
    if not 0 <= int(lesion_class) < labels.size:
        raise ValueError("lesion_class is outside template_labels")

    best_score = np.zeros(ordinary.shape[:3], dtype=np.float32)
    best_class = np.full(ordinary.shape[:3], -1, dtype=np.intp)
    for class_index in range(ordinary.shape[-1]):
        score = ordinary[..., class_index]
        selected = score > best_score
        best_class[selected] = class_index
        best_score[selected] = score[selected]

    result = np.zeros(ordinary.shape[:3], dtype=np.int8)
    ordinary_wins = best_class >= 0
    result[ordinary_wins] = labels.astype(np.int8, copy=False)[best_class[ordinary_wins]]
    lesion_wins = lesion > best_score
    result[lesion_wins] = np.int8(labels[int(lesion_class)])
    return result
