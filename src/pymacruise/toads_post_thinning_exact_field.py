"""One selected next-loop TOADS degree-one exact field update."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_inhomogeneity_exact_fit import (
    ToadsExactChebyshevFieldFit,
    fit_toads_exact_chebyshev_field_degree1,
)
from .toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from .toads_post_registration_next_thinning import ToadsPostRegistrationNextThinning


@dataclass(frozen=True)
class ToadsPostThinningExactField:
    """The next thinning state and the field-refit stage selected by it."""

    next_thinning: ToadsPostRegistrationNextThinning
    degree: int
    fit: ToadsExactChebyshevFieldFit
    stage_order: tuple[str, ...]


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
    fit = fit_toads_exact_chebyshev_field_degree1(
        image_values,
        branch.handoff.membership_pass.memberships,
        next_thinning.thinning.segmentation,
        centroids,
        field,
    )
    return ToadsPostThinningExactField(
        next_thinning=next_thinning,
        degree=degree,
        fit=fit,
        stage_order=(
            *next_thinning.stage_order,
            "next_iteration.computeExact3DPolynomialField.degree1",
        ),
    )
