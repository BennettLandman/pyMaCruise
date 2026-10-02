"""Degree-gated exact field updates after a selected TOADS thinning pass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_inhomogeneity_exact_fit import (
    ToadsExactChebyshevFieldFit,
    fit_toads_exact_chebyshev_field,
)
from .toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from .toads_post_registration_next_thinning import ToadsPostRegistrationNextThinning


@dataclass(frozen=True)
class ToadsPostThinningExactField:
    """The next thinning state and the field-refit stage selected by it."""

    next_thinning: ToadsPostRegistrationNextThinning
    previous_degree: int
    degree: int
    fit: ToadsExactChebyshevFieldFit | None
    stage_order: tuple[str, ...]


def toads_post_thinning_exact_field_update(
    next_thinning: ToadsPostRegistrationNextThinning,
    image_values: Any,
    centroids: Any,
    field: Any,
    *,
    current_degree: int,
    polynomial_degree: int = 3,
    max_distance: float = 0.0001,
    low_energy_ratio: float | None = None,
) -> ToadsPostThinningExactField:
    """Apply the source one-degree ramp and exact field update/no-op.

    `AlgorithmMultipleToads` starts with degree zero and increments by one
    when its previous energy ratio is strictly below `sqrt(maxDistance)`.
    `computeCorrectionField(0)` returns without fitting; every positive
    degree refits the attached field, including an unchanged maximum degree,
    using the exact-membership Chebyshev method.
    """
    if not isinstance(next_thinning, ToadsPostRegistrationNextThinning):
        raise TypeError("next_thinning must be a selected post-registration thinning result")
    branch = next_thinning.branch
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("field update requires a continuing registered branch")
    if isinstance(max_distance, bool) or not np.isfinite(np.float32(max_distance)) or np.float32(max_distance) < 0:
        raise ValueError("max_distance must be finite and nonnegative")
    if low_energy_ratio is not None and (
        isinstance(low_energy_ratio, bool)
        or not np.isfinite(np.float32(low_energy_ratio))
        or np.float32(low_energy_ratio) < 0
    ):
        raise ValueError("low_energy_ratio must be finite and nonnegative")
    previous_degree = current_degree
    degree = toads_next_inhomogeneity_degree(
        current_degree,
        polynomial_degree,
        branch.control.energy_ratio,
        np.float32(low_energy_ratio)
        if low_energy_ratio is not None
        else np.sqrt(np.float32(max_distance)),
    )
    if degree == 0:
        return ToadsPostThinningExactField(
            next_thinning=next_thinning,
            previous_degree=previous_degree,
            degree=degree,
            fit=None,
            stage_order=(
                *next_thinning.stage_order,
                "next_iteration.computeExact3DPolynomialField.degree0.noop",
            ),
        )
    fit = fit_toads_exact_chebyshev_field(
        image_values,
        branch.handoff.membership_pass.memberships,
        next_thinning.thinning.segmentation,
        centroids,
        field,
        degree=degree,
    )
    return ToadsPostThinningExactField(
        next_thinning=next_thinning,
        previous_degree=previous_degree,
        degree=degree,
        fit=fit,
        stage_order=(
            *next_thinning.stage_order,
            f"next_iteration.computeExact3DPolynomialField.degree{degree}",
        ),
    )


def toads_post_thinning_exact_field_degree1(
    next_thinning: ToadsPostRegistrationNextThinning,
    image_values: Any,
    centroids: Any,
    field: Any,
) -> ToadsPostThinningExactField:
    """Advance a continuing selected branch through its degree-one exact fit.

    The source starts each short-constructor run at degree zero, raises it by
    one when the just-recorded energy ratio is below ``sqrt(0.0001)``, then
    performs one global field fit per channel. This bounded adapter handles
    exactly the selected single-channel degree-one case; other degrees and
    multichannel scheduling are deliberately not inferred here.
    """
    if not isinstance(next_thinning, ToadsPostRegistrationNextThinning):
        raise TypeError("next_thinning must be a selected post-registration thinning result")
    branch = next_thinning.branch
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("degree-one field update requires a continuing registered branch")
    degree = toads_next_inhomogeneity_degree(
        0,
        3,
        branch.control.energy_ratio,
        np.sqrt(np.float32(0.0001)),
    )
    if degree != 1:
        raise ValueError("the selected branch energy ratio does not schedule degree one")
    result = toads_post_thinning_exact_field_update(
        next_thinning,
        image_values,
        centroids,
        field,
        current_degree=0,
    )
    if result.degree != 1 or result.fit is None:
        raise ValueError("the selected branch energy ratio does not schedule degree one")
    return result
