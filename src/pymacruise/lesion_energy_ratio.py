"""The selected LesionTOADS wrapper's float32 energy-ratio expression.

This is only the scalar state update after ``computeMemberships()``. It does
not decide convergence or perform registration/growth. Arithmetic and
``Numerics.abs(float)`` follow the selected Java expression, including the
source method's sign flip for NaNs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LesionEnergyRatio:
    """Float32 operands and intermediate values for one wrapper update."""

    previous_energy: np.float32
    current_energy: np.float32
    difference: np.float32
    numerator: np.float32
    denominator: np.float32
    quotient: np.float32
    energy_ratio: np.float32


def lesion_toads_energy_ratio(
    previous_energy: float,
    current_energy: float,
) -> LesionEnergyRatio:
    """Evaluate ``Numerics.abs(2.0f*(previous-current)/(previous+current))``.

    Explicit float32 materialization preserves Java's left-to-right float
    operations. The final absolute value mirrors VandyCRUISE
    ``Numerics.abs(float)`` (``a > 0 ? a : -a``), which differs from
    ``Math.abs`` for the sign bit of a NaN. IEEE exceptional values are
    preserved because the selected membership pass can return NaN.
    """

    previous = np.float32(previous_energy)
    current = np.float32(current_energy)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        difference = np.float32(previous - current)
        numerator = np.float32(np.float32(2.0) * difference)
        denominator = np.float32(previous + current)
        quotient = np.float32(numerator / denominator)
        ratio = np.float32(quotient if quotient > np.float32(0.0) else -quotient)
    return LesionEnergyRatio(
        previous_energy=previous,
        current_energy=current,
        difference=difference,
        numerator=numerator,
        denominator=denominator,
        quotient=quotient,
        energy_ratio=ratio,
    )
