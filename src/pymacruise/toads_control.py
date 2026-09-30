"""Selected outer-loop control decisions for AlgorithmMultipleToads."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


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
