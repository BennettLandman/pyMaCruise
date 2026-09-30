"""Selected base-TOADS thalamus-membership output adapter.

This module ports only ``MultipleToadSegmentation.generateThalamusOutput``.
It does not evolve thalamus labels or perform MGDM segmentation. Volumes use
the selected Java ``[x,y,z]`` order, with class/candidate axes last; image
orientation and persistence are handled by the caller.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def generate_thalamus_output(
    memberships: ArrayLike,
    atlas_names: tuple[str, ...] | list[str],
    *,
    use_approx: bool,
    best: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Return the one-channel output from selected Java's thalamus adapter.

    Exact mode copies the membership volume for every atlas entry named exactly
    ``"Thalamus"``. Thus a later duplicate name overwrites the earlier one.
    Approximate mode processes atlas classes in the same order, and candidate
    slots in ascending order within each class; every matching slot overwrites
    the prior output, so the last matching slot of the last duplicate class
    wins. Missing and case/whitespace-near-match names leave the zero buffer.

    The returned array has shape ``(1, x, y, z)`` and dtype float32, matching
    Java's ``float[1][nx][ny][nz]`` result. Use
    :func:`pymacruise.toads_crop.uncrop_toads_buffer` for wrapper uncropping.
    """
    mem = np.asarray(memberships)
    names = tuple(atlas_names)
    if mem.ndim != 4 or any(size < 1 for size in mem.shape[:3]):
        raise ValueError("memberships must have shape [x,y,z,class]")
    if not (
        np.issubdtype(mem.dtype, np.integer)
        or np.issubdtype(mem.dtype, np.floating)
    ):
        raise TypeError("memberships must be real numeric values")
    if len(names) == 0 or mem.shape[3] != len(names):
        raise ValueError("membership class axis must match a nonempty atlas name list")
    if any(not isinstance(name, str) for name in names):
        raise TypeError("atlas_names must contain strings")
    if not isinstance(use_approx, (bool, np.bool_)):
        raise TypeError("use_approx must be boolean")

    candidate_classes: NDArray[np.integer] | None = None
    if use_approx:
        if best is None:
            raise ValueError("best candidate class indices are required in approximate mode")
        candidate_classes = np.asarray(best)
        if candidate_classes.ndim != 4 or candidate_classes.shape[:3] != mem.shape[:3]:
            raise ValueError("best must have shape [x,y,z,candidate]")
        if candidate_classes.shape[3] < 1 or not np.issubdtype(candidate_classes.dtype, np.integer):
            raise ValueError("best must contain one or more integer candidate class indices")
        if candidate_classes.shape[3] > mem.shape[3]:
            raise ValueError("best candidate slots exceed the membership candidate axis")
        if np.any(candidate_classes < 0) or np.any(candidate_classes >= len(names)):
            raise ValueError("best contains an atlas class index outside atlas_names")

    memberships32 = np.asarray(mem, dtype=np.float32)
    output = np.zeros((1, *mem.shape[:3]), dtype=np.float32)
    for class_index, name in enumerate(names):
        if name != "Thalamus":
            continue
        if candidate_classes is None:
            output[0, ...] = memberships32[..., class_index]
        else:
            for slot in range(candidate_classes.shape[3]):
                selected = candidate_classes[..., slot] == class_index
                output[0, selected] = memberships32[..., slot][selected]
    return output
