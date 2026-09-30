"""Selected VandyCRUISE TOADS hard-classification output helper.

This reproduces ``MultipleToadSegmentation.generateClassification`` for its
two class-selection routes. It returns signed Java-byte labels in the same
``[x, y, z]`` geometry as the class arrays. Crop restoration and serialization
remain separate concerns; callers can pass the result to the already validated
``uncrop_toads_buffer`` adapter when they need Java's flattened buffer layout.
"""

from __future__ import annotations

import numpy as np


def generate_toads_classification(
    template_labels: np.ndarray,
    memberships: np.ndarray,
    *,
    use_approx: bool = False,
    best_class_indices: np.ndarray | None = None,
) -> np.ndarray:
    """Map winning class indices to the selected source's template byte labels.

    ``memberships`` has shape ``[x, y, z, class]``. In exact mode, the source
    scans classes in ascending order with ``max=0`` and updates only on strict
    ``>``; ties therefore retain the first winner and all nonpositive/NaN
    scores leave the output at byte zero. Positive infinity compares greater
    than the initial maximum and may win. In approximate mode, the source uses
    only candidate ``best[..., 0]`` and ignores the membership values.

    Java ``byte`` is signed. The returned NumPy ``int8`` values preserve the
    signed representation of each template label; they are labels, not class
    indices. No additional saturation or unsigned reinterpretation is applied.
    """
    template = np.asarray(template_labels)
    mem = np.asarray(memberships)
    if template.ndim != 1 or not np.issubdtype(template.dtype, np.integer):
        raise ValueError("template_labels must be a one-dimensional integer byte-label array")
    if np.any(template < -128) or np.any(template > 127):
        raise ValueError("template_labels must contain signed Java byte values in [-128, 127]")
    if mem.ndim != 4 or not np.issubdtype(mem.dtype, np.number) or np.iscomplexobj(mem):
        raise ValueError("memberships must be a real numeric [x, y, z, class] array")
    # The selected Java field is float[][][][]; compare after Java-float
    # precision conversion even if a caller supplied a wider NumPy dtype.
    mem = mem.astype(np.float32, copy=False)
    if mem.shape[-1] != template.size or template.size == 0:
        raise ValueError("the membership class axis must match nonempty template_labels")
    shape = mem.shape[:3]

    if use_approx:
        if best_class_indices is None:
            raise ValueError("best_class_indices are required when use_approx is true")
        best = np.asarray(best_class_indices)
        if best.shape != shape or not np.issubdtype(best.dtype, np.integer):
            raise ValueError("best_class_indices must be an integer [x, y, z] array")
        if np.any(best < 0) or np.any(best >= template.size):
            raise ValueError("best_class_indices contains an index outside template_labels")
        return template.astype(np.int8, copy=False)[best]

    # Keep the source's strict scalar comparison semantics. In particular,
    # np.argmax would incorrectly select a zero/negative/NaN-only class.
    max_score = np.zeros(shape, dtype=mem.dtype)
    winner = np.full(shape, -1, dtype=np.intp)
    for class_index in range(template.size):
        score = mem[..., class_index]
        wins = score > max_score
        winner[wins] = class_index
        max_score[wins] = score[wins]

    result = np.zeros(shape, dtype=np.int8)
    selected = winner >= 0
    result[selected] = template.astype(np.int8, copy=False)[winner[selected]]
    return result
