"""Run selected ``registerShapes`` after a 13-class GM membership update."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_next_iteration_membership import LesionNextIterationMembership
from .lesion_registration_handoff import (
    LesionRegistrationHandoff,
    register_lesion_membership_handoff,
)


@dataclass(frozen=True)
class LesionNextIterationRegistration:
    """Membership update, energy ratio, and selected atlas-registration handoff."""

    membership_update: LesionNextIterationMembership
    registration: LesionRegistrationHandoff
    stage_order: tuple[str, ...]


def register_after_next_gm_membership_update(
    membership_update: LesionNextIterationMembership,
    shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    align_shapes: bool = True,
    init_align_iter: int = 10,
) -> LesionNextIterationRegistration:
    """Call selected ``registerShapes()`` immediately after membership energy.

    The selected ``AlgorithmLesionToads`` wrapper performs atlas registration
    before evaluating the energy-ratio convergence check. Therefore the ratio
    is retained as state and validated against this membership objective, but
    it does not gate this registration call.
    """
    if not isinstance(membership_update, LesionNextIterationMembership):
        raise TypeError("membership_update must be a LesionNextIterationMembership")
    ratio = membership_update.energy_ratio
    objective = np.float32(membership_update.membership.objective)
    current = np.float32(ratio.current_energy)
    if objective.view(np.uint32) != current.view(np.uint32):
        raise ValueError("energy ratio current_energy must match membership objective bits")
    registration = register_lesion_membership_handoff(
        membership_update.membership,
        shape_priors,
        initial_transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    return LesionNextIterationRegistration(
        membership_update,
        registration,
        membership_update.stage_order + ("registerShapes",),
    )
