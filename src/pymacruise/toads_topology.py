"""Selected VandyCRUISE multi-label 18/6 homeomorphism gate."""

from __future__ import annotations

import numpy as np

from pymacruise.toads_critical_lut import critical186_is_regular


def toads_homeomorphic_label_change(
    neighborhood: np.ndarray,
    candidate_label: int,
) -> bool:
    """Return whether one center-label change passes TOADS' 18/6 LUT checks.

    The input is a complete source x/y/z ordered 3×3×3 label neighborhood;
    the center currently has its old label and ``candidate_label`` is the
    proposed replacement. This ports the selected
    `MultipleToadSegmentation.isHomeomorphicSegmentation` check with its
    constructor setting ``checkComposed == false``: test the new object, old
    object, each old/new neighboring-label pair, and each old/new triplet.

    It is a local admissibility predicate only. It does not evolve memberships,
    choose labels, check out-of-volume borders, or reproduce segmentation
    update ordering. The caller must provide an interior neighborhood.
    """
    labels = np.asarray(neighborhood)
    if labels.shape != (3, 3, 3):
        raise ValueError("TOADS topology neighborhood must have shape (3, 3, 3)")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("TOADS topology labels must be integers")
    if np.any(labels < -32768) or np.any(labels > 32767):
        raise ValueError("TOADS topology labels must fit Java's signed short range")
    if not isinstance(candidate_label, (int, np.integer)) or not -32768 <= candidate_label <= 32767:
        raise ValueError("candidate_label must fit Java's signed short range")

    old_label = int(labels[1, 1, 1])
    new_label = int(candidate_label)
    if old_label == new_label:
        return True

    # MultipleToadSegmentation builds a separate boolean object image for each
    # relation, setting the changed center to its new/old state as appropriate.
    new_object = labels == new_label
    new_object[1, 1, 1] = True
    if not _regular(new_object):
        return False

    old_object = labels == old_label
    old_object[1, 1, 1] = False
    if not _regular(old_object):
        return False

    # The Java loop gathers distinct neighboring labels in x/y/z scan order,
    # excluding the proposed new label and the old center label.
    related_labels: list[int] = []
    for x in range(3):
        for y in range(3):
            for z in range(3):
                if (x, y, z) == (1, 1, 1):
                    continue
                label = int(labels[x, y, z])
                if label in (new_label, old_label) or label in related_labels:
                    continue
                related_labels.append(label)

    for relation in related_labels:
        old_pair = (labels == old_label) | (labels == relation)
        old_pair[1, 1, 1] = False
        if not _regular(old_pair):
            return False

        new_pair = (labels == new_label) | (labels == relation)
        new_pair[1, 1, 1] = True
        if not _regular(new_pair):
            return False

    for first_index, first in enumerate(related_labels):
        for second in related_labels[first_index + 1 :]:
            old_triplet = (labels == old_label) | (labels == first) | (labels == second)
            old_triplet[1, 1, 1] = False
            if not _regular(old_triplet):
                return False

            new_triplet = (labels == new_label) | (labels == first) | (labels == second)
            new_triplet[1, 1, 1] = True
            if not _regular(new_triplet):
                return False

    return True


def _regular(neighborhood: np.ndarray) -> bool:
    """Apply the selected source LUT lookup to one binary relation."""
    return critical186_is_regular(neighborhood)
