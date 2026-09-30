"""Selected LesionTOADS membership-derived lesion classification map.

This isolates ``LesionToadSegmentation.updateLesionClassifcationformMemberships``
from the selected VandyCRUISE source. It returns the short-valued
``lesionClassification_mem`` map and its source-counted volume; it does not
claim to be an active LesionTOADS workflow stage because the supplied source
tree contains no ``PlugInLesionToads`` caller for this method.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LesionMembershipClassification:
    """Binary Java-short-compatible map and source lesion-volume count."""

    classification_map: NDArray[np.int16]
    lesion_volume: int
    gm_class_index: int


def lesion_toads_membership_classification(
    classification: ArrayLike,
    template_labels: ArrayLike,
    atlas_names: Sequence[str],
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    lesion_class: int,
    *,
    has_lesion: bool = True,
) -> LesionMembershipClassification:
    """Build the membership-derived map using the selected Java comparisons.

    Atlas class names are searched case-insensitively for
    ``"Cerebellum-GM"``; the last match wins, and the Java default index 0 is
    retained when there is no match. A voxel is marked one exactly when
    ``has_lesion`` is true, its signed-byte classification equals the Java
    byte label at ``template_labels[lesion_class]``, and its float32 ordinary
    lesion-class membership is strictly less than the float32 sum of lesion
    membership and the selected GM membership. The sum is rounded to float32
    before the strict comparison, as Java float addition does.

    Arrays use the package's RAS-oriented ``(x, y, z)`` voxel convention;
    memberships append the class axis. The returned map has Java-short
    compatible ``int16`` storage and contains only zero and one.
    """
    labels = np.asarray(classification)
    template = np.asarray(template_labels)
    class_memberships = np.asarray(memberships, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)

    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("classification must be a nonempty 3-D (x, y, z) array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise TypeError("classification must contain Java byte labels")
    if np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("classification values must fit Java's signed byte range")
    if template.ndim != 1 or not np.issubdtype(template.dtype, np.integer):
        raise TypeError("template_labels must be a one-dimensional integer array")
    if np.any(template < -128) or np.any(template > 127):
        raise ValueError("template_labels values must fit Java's signed byte range")
    if class_memberships.ndim != 4 or class_memberships.shape[:3] != labels.shape:
        raise ValueError("memberships must match classification with a class axis")
    classes = class_memberships.shape[-1]
    if isinstance(lesion_class, (bool, np.bool_)) or not isinstance(
        lesion_class, (int, np.integer)
    ):
        raise ValueError("lesion_class must be an integer class index")
    if not 0 <= int(lesion_class) < classes:
        raise ValueError("lesion_class is outside the membership class axis")
    if classes != len(template):
        raise ValueError("template_labels and memberships must have the same class count")
    if lesion.shape != labels.shape:
        raise ValueError("lesion_membership must match classification")
    if not np.isfinite(class_memberships).all() or not np.isfinite(lesion).all():
        raise ValueError("memberships must contain only finite values")
    if not isinstance(has_lesion, (bool, np.bool_)):
        raise ValueError("has_lesion must be boolean")

    names = tuple(atlas_names)
    if len(names) != classes or any(not isinstance(name, str) for name in names):
        raise ValueError("atlas_names must contain one string per membership class")

    gm_class = 0
    for index, name in enumerate(names):
        if name.casefold() == "cerebellum-gm".casefold():
            gm_class = index

    result = np.zeros(labels.shape, dtype=np.int16)
    if has_lesion:
        target = np.int8(template[int(lesion_class)])
        ordinary = class_memberships[..., int(lesion_class)]
        # Explicit float32 output preserves Java's float + float rounding.
        lesion_plus_gm = np.add(lesion, class_memberships[..., gm_class], dtype=np.float32)
        selected = (labels.astype(np.int8, copy=False) == target) & (ordinary < lesion_plus_gm)
        result[selected] = 1

    return LesionMembershipClassification(
        classification_map=result,
        lesion_volume=int(np.count_nonzero(result)),
        gm_class_index=gm_class,
    )
