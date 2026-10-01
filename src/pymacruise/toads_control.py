"""Selected outer-loop control decisions for AlgorithmMultipleToads.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
TOADS as an iterative statistical and topology-preserving segmentation
framework. This module implements only the selected wrapper's stop, growth,
classification, and relation-refresh decisions after one membership-energy
update. Thresholds and stop ordering are source-specific; this control slice
does not implement or validate paper-level convergence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class ToadsIterationControl:
    """Control choices after one membership-energy update."""

    energy_ratio: np.float32
    stop: bool
    run_growing: bool
    update_classification: bool
    refresh_relations: bool


def toads_iteration_control(
    previous_energy: float,
    current_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    *,
    do_growing: bool,
) -> ToadsIterationControl:
    """Return source-ordered stop, growing, classification, and refresh choices.

    This ports only the `AlgorithmMultipleToads` decisions after
    `computeMemberships()`. The energy change uses Java float32 subtraction,
    multiplication, addition, division, and absolute value in the source
    expression order. The selected source still runs growing on the stopping
    iteration when ``do_growing`` is enabled.
    """
    if isinstance(iteration, (bool, np.bool_)) or not isinstance(iteration, (int, np.integer)) or iteration < 1:
        raise ValueError("iteration must be an integer greater than or equal to one")
    if (
        isinstance(maximum_iterations, (bool, np.bool_))
        or not isinstance(maximum_iterations, (int, np.integer))
        or maximum_iterations < 1
    ):
        raise ValueError("maximum_iterations must be an integer greater than or equal to one")
    if not isinstance(do_growing, (bool, np.bool_)):
        raise ValueError("do_growing must be boolean")
    previous = np.float32(previous_energy)
    current = np.float32(current_energy)
    threshold = np.float32(convergence_threshold)
    if not np.isfinite(previous) or not np.isfinite(current):
        raise ValueError("membership energies must be finite")
    if not np.isfinite(threshold) or threshold < np.float32(0.0):
        raise ValueError("convergence_threshold must be finite and nonnegative")

    difference = np.float32(previous - current)
    numerator = np.float32(np.float32(2.0) * difference)
    denominator = np.float32(previous + current)
    ratio = np.float32(np.abs(np.float32(numerator / denominator)))
    stop = bool(ratio < threshold)
    if iteration >= maximum_iterations:
        stop = True
    run_growing = bool((not stop) or do_growing)
    return ToadsIterationControl(
        energy_ratio=ratio,
        stop=stop,
        run_growing=run_growing,
        update_classification=stop,
        refresh_relations=not stop,
    )


def toads_classification_from_labels(
    segmentation: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
) -> tuple[NDArray[np.int8], int]:
    """Map segmentation class indices to source signed-byte output labels.

    This follows ``MultipleToadSegmentation.updateClassificationFromLabels``:
    visit every voxel, assign ``templateLabel[segmentation]``, and count
    changed output bytes. Membership winner selection is a separate source
    method and is intentionally not performed here.
    """
    labels = np.asarray(segmentation)
    template = np.asarray(template_labels)
    previous = np.asarray(previous_classification)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if template.ndim != 1 or template.size == 0:
        raise ValueError("template_labels must be a nonempty one-dimensional array")
    if not np.issubdtype(template.dtype, np.integer) or np.issubdtype(template.dtype, np.bool_):
        raise ValueError("template_labels must be integers")
    if previous.shape != labels.shape:
        raise ValueError("previous_classification must match segmentation")
    if not np.issubdtype(previous.dtype, np.integer) or np.issubdtype(previous.dtype, np.bool_):
        raise ValueError("previous_classification must contain integer byte labels")
    if np.any(labels < 0) or np.any(labels >= template.size):
        raise ValueError("segmentation contains an out-of-range template-label index")
    byte_limits = np.iinfo(np.int8)
    if np.any(template < byte_limits.min) or np.any(template > byte_limits.max):
        raise ValueError("template_labels must fit in Java signed-byte range")
    if np.any(previous < byte_limits.min) or np.any(previous > byte_limits.max):
        raise ValueError("previous_classification must fit in Java signed-byte range")

    output = template[labels].astype(np.int8, copy=False)
    changed = int(np.count_nonzero(previous != output))
    return output, changed
