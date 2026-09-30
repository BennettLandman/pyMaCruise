"""Composed post-membership stopping transition from AlgorithmMultipleToads."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_control import (
    ToadsIterationControl,
    toads_classification_from_labels,
    toads_iteration_control,
)
from pymacruise.toads_thinning import ToadsGrowingPass, toads_propagate_growing


@dataclass(frozen=True)
class ToadsStoppingTransition:
    """Selected final iteration's growth and classification result."""

    control: ToadsIterationControl
    growing: ToadsGrowingPass | None
    classification: NDArray[np.int8]
    changed_classification_voxels: int


def toads_stopping_iteration(
    previous_energy: float,
    current_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    *,
    do_growing: bool,
    available: ArrayLike | None = None,
) -> ToadsStoppingTransition:
    """Apply the selected stop branch after membership/energy calculation.

    This composes the source's energy/iteration stop decision, optional final
    ``propagateGrowing()``, and ``updateClassificationFromLabels()``. The
    selected main loop does not refresh relation maps after it stops. Thinning,
    field correction, centroid/variance estimation, membership updates, and
    atlas registration remain separate preceding stages.
    """
    control = toads_iteration_control(
        previous_energy,
        current_energy,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    if not control.stop:
        raise ValueError("toads_stopping_iteration requires a selected stop condition")

    labels = np.asarray(segmentation)
    skeleton_values = np.asarray(skeleton)
    order_values = np.asarray(ordering, dtype=np.float32)
    volumes = np.asarray(class_volumes, dtype=np.float32)
    growing = None
    if control.run_growing:
        growing = toads_propagate_growing(
            labels,
            skeleton_values,
            memberships,
            order_values,
            volumes,
            initial_available=available,
        )
        labels = growing.segmentation

    classification, changed = toads_classification_from_labels(
        labels, template_labels, previous_classification
    )
    return ToadsStoppingTransition(control, growing, classification, changed)
