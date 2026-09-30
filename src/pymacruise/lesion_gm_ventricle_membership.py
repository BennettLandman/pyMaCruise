"""Selected GM-class ventricle attenuation in LesionTOADS memberships.

This is a bounded formula excerpt from ``LesionToadSegmentation.computeMemberships``.
It assumes the caller has already calculated each class's scalar ``num`` term;
it applies the source's exact class-name gate, DIRECT/PRIOR expression order,
denominator accumulation, and final ordinary-class normalization for one voxel.
It is not a complete membership estimator.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray


_VENTRICLE_ATTENUATED_GM_NAMES = frozenset({"Cerebrum-GM", "GM"})


def lesion_toads_gm_ventricle_memberships(
    class_names: Sequence[str],
    membership_terms: np.ndarray,
    shape_values: np.ndarray,
    *,
    relation_type: str,
    has_lesion: bool,
    ventricle_distance: int,
    max_ventricle_distance: int,
) -> NDArray[np.float32]:
    """Apply the selected one-voxel ordinary-class membership update.

    ``membership_terms`` are the already-computed Java ``num`` values. Classes
    are visited in the supplied order. Under ``has_lesion``, only the exact
    names ``Cerebrum-GM`` and ``GM`` are multiplied by signed-byte
    ``ventricle_distance / max_ventricle_distance`` when the maximum is
    nonzero and the distance is not ``-1``. Otherwise they follow the ordinary
    DIRECT or PRIOR assignment. PRIOR applies ``min(1, 2*shape)`` before the
    attenuation, preserving Java's float32 expression order. The ordinary
    denominator is accumulated in class order and positive denominators are
    applied in a second pass, matching the selected excerpt.

    The inputs describe one voxel in the source's RAS-oriented class ordering;
    this scalar formula does not perform image reorientation.
    """
    names = tuple(class_names)
    terms = np.asarray(membership_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if not names or terms.ndim != 1 or shapes.ndim != 1:
        raise ValueError("class names, membership_terms, and shape_values must be nonempty vectors")
    if terms.shape != (len(names),) or shapes.shape != (len(names),):
        raise ValueError("membership_terms and shape_values must have one value per class")
    if any(not isinstance(name, str) for name in names):
        raise TypeError("class_names must contain only strings")
    if relation_type not in {"DIRECT", "PRIOR"}:
        raise ValueError("relation_type must be DIRECT or PRIOR")
    if not np.isfinite(terms).all() or not np.isfinite(shapes).all():
        raise ValueError("membership_terms and shape_values must be finite")
    if isinstance(ventricle_distance, bool) or not isinstance(ventricle_distance, (int, np.integer)):
        raise TypeError("ventricle_distance must be a signed integer distance")
    if isinstance(max_ventricle_distance, bool) or not isinstance(max_ventricle_distance, (int, np.integer)):
        raise TypeError("max_ventricle_distance must be an integer")

    distance = np.float32(ventricle_distance)
    maximum = np.float32(max_ventricle_distance)
    attenuate = bool(has_lesion and maximum != np.float32(0.0) and distance != np.float32(-1.0))
    result = np.empty(len(names), dtype=np.float32)
    denominator = np.float32(0.0)

    for index, name in enumerate(names):
        num = terms[index]
        if relation_type == "PRIOR":
            doubled_shape = np.float32(np.float32(2.0) * shapes[index])
            shape_factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
            value = np.float32(shape_factor * num)
        else:
            value = np.float32(num)

        if attenuate and name in _VENTRICLE_ATTENUATED_GM_NAMES:
            if relation_type == "PRIOR":
                value = np.float32(value * distance)
                value = np.float32(value / maximum)
            else:
                value = np.float32(num * distance)
                value = np.float32(value / maximum)

        result[index] = value
        denominator = np.float32(denominator + value)

    if denominator > np.float32(0.0):
        for index in range(len(names)):
            result[index] = np.float32(result[index] / denominator)
    return result
