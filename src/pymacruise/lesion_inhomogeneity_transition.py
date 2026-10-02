"""Selected LesionTOADS post-thinning inhomogeneity degree/refit transition.

``AlgorithmLesionToads.runAlgorithm`` raises the correction degree by at most
one when the previous energy ratio is strictly below ``sqrt(maxDistance)``;
it then calls ``computeCorrectionField(deg)`` once per channel. The refit is
injected so the selected ``InhomogeneityCorrection2`` kernel can be reused
without duplicating its exact/approximate numerical implementation here.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


ChannelRefit = Callable[[int, NDArray[np.float32], Any], Any]


@dataclass(frozen=True)
class LesionToadsInhomogeneityTransition:
    """Degree state, attached fields, and per-channel fit diagnostics."""

    previous_degree: int
    degree: int
    low_energy_ratio: np.float32
    fields: tuple[NDArray[np.float32], ...]
    fits: tuple[Any, ...]
    stage_order: tuple[str, ...]


def update_lesion_toads_inhomogeneity_fields(
    attached_fields: Sequence[ArrayLike],
    refit_channel: ChannelRefit,
    *,
    degree: int,
    polynomial_degree: int,
    previous_energy_ratio: float,
    max_distance: float,
) -> LesionToadsInhomogeneityTransition:
    """Apply the source degree gate and invoke each channel's correction.

    ``maxDistance`` is converted to Java float before ``Math.sqrt`` and its
    result is rounded back to float32, matching ``lowEratio`` in the wrapper.
    NaN energy ratios retain Java comparison behavior (the strict ``<`` is
    false), while infinities are rejected as unsupported inputs. Each
    callback receives its channel index, the same attached float32 field
    object supplied by the caller, and the selected degree. Degree zero still
    calls every channel, matching ``computeCorrectionField(0)`` (a no-op in
    the selected Java correction implementation).
    """
    for name, value in (("degree", degree), ("polynomial_degree", polynomial_degree)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer")
        if int(value) < 0:
            raise ValueError(f"{name} must be nonnegative")
    if int(degree) > int(polynomial_degree):
        raise ValueError("degree cannot exceed polynomial_degree")
    if not callable(refit_channel):
        raise TypeError("refit_channel must be callable")
    max_distance_value = np.float32(max_distance)
    ratio = np.float32(previous_energy_ratio)
    if not np.isfinite(max_distance_value) or max_distance_value < np.float32(0.0):
        raise ValueError("max_distance must be finite and nonnegative")
    if np.isinf(ratio):
        raise ValueError("previous_energy_ratio cannot be infinite")

    # Java evaluates (float)Math.sqrt(float maxDistance); preserve that
    # double-then-float sequence instead of relying on a float32 ufunc.
    low_ratio = np.float32(math.sqrt(float(max_distance_value)))
    selected_degree = int(degree)
    if selected_degree < int(polynomial_degree) and ratio < low_ratio:
        selected_degree += 1

    try:
        fields = tuple(np.asarray(field) for field in attached_fields)
    except TypeError as exc:
        raise ValueError("attached_fields must be a sequence of channel fields") from exc
    if not fields:
        raise ValueError("at least one attached correction field is required")
    if any(field.ndim != 3 or field.dtype != np.float32 for field in fields):
        raise ValueError("each attached correction field must be a 3-D float32 array")
    if any(field.shape != fields[0].shape for field in fields[1:]):
        raise ValueError("attached correction fields must share the same grid")

    fits: list[Any] = []
    stages = ["previousEnergyRatio", "degreeGate"]
    for channel, field in enumerate(fields):
        fit = refit_channel(channel, field, selected_degree)
        fits.append(fit)
        stages.append(f"computeCorrectionField.channel{channel}.degree{selected_degree}")
    return LesionToadsInhomogeneityTransition(
        previous_degree=int(degree),
        degree=selected_degree,
        low_energy_ratio=low_ratio,
        fields=fields,
        fits=tuple(fits),
        stage_order=tuple(stages),
    )
