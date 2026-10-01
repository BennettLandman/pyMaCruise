"""Call the selected LesionTOADS atlas registration after membership updates.

This source-ordered continuation is ``computeMemberships`` → energy ratio →
``registerShapes`` from the VandyCRUISE linear-fast-marching wrapper. It stops
at the registration handoff and does not make a convergence or output claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_post_thinning_membership import LesionPostThinningMembershipUpdate
from .lesion_energy_ratio import LesionEnergyRatio
from .lesion_registration_handoff import (
    LesionRegistrationHandoff,
    register_lesion_membership_handoff,
)


@dataclass(frozen=True)
class LesionPostThinningRegistration:
    """Membership, wrapper ratio, and the immediately following registration."""

    membership_update: LesionPostThinningMembershipUpdate
    energy_ratio: LesionEnergyRatio
    registration: LesionRegistrationHandoff
    stage_order: tuple[str, ...]


def register_after_post_thinning_membership(
    membership_update: LesionPostThinningMembershipUpdate,
    energy_ratio: LesionEnergyRatio,
    shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    align_shapes: bool = True,
    init_align_iter: int = 10,
) -> LesionPostThinningRegistration:
    """Run selected ``registerShapes()`` after P39am/P39an state.

    Membership and the energy-ratio transition are supplied completed stages.
    The ratio is deliberately not used as a gate: the selected wrapper invokes
    ``atlas.registerShapes()`` before its convergence check. Consequently a
    source-produced NaN ratio does not suppress this registration call, while
    the membership arrays themselves still must be finite for optimization.
    """
    if not isinstance(membership_update, LesionPostThinningMembershipUpdate):
        raise TypeError("membership_update must be a LesionPostThinningMembershipUpdate")
    if not isinstance(energy_ratio, LesionEnergyRatio):
        raise TypeError("energy_ratio must be a LesionEnergyRatio")
    current = np.float32(membership_update.membership.objective)
    supplied = np.float32(energy_ratio.current_energy)
    # Compare raw bits so the selected NaN payload/sign from Java is accepted
    # only when it corresponds to the objective of the supplied membership.
    if current.view(np.uint32) != supplied.view(np.uint32):
        raise ValueError("energy_ratio current_energy must match membership objective bits")

    registration = register_lesion_membership_handoff(
        membership_update.membership,
        shape_priors,
        initial_transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    return LesionPostThinningRegistration(
        membership_update,
        energy_ratio,
        registration,
        membership_update.stage_order + ("energyRatio", "registerShapes"),
    )
