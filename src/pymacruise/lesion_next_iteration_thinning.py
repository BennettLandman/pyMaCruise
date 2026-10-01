"""Handoff a continuing LesionTOADS branch to its next thinning pass.

This isolates the selected default-wrapper loop restart:
``propagateGrowing`` has produced a frontier, and the next iteration starts
with ``propagateThinning``. Statistics, memberships, registration, and the
following energy pass are outside this adapter. The input profile is explicit
so the operation can accept the selected 13-class MaCRUISEplus state as well
as smaller source-harness profiles.
"""

from __future__ import annotations

from typing import Sequence

from numpy.typing import ArrayLike

from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_post_registration import LesionToadsPostRegistrationBranch


def run_lesion_toads_next_iteration_thinning(
    previous_branch: LesionToadsPostRegistrationBranch,
    transformed_shape: ArrayLike,
    *,
    object_types: Sequence[str],
    lesion_class: int,
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
) -> LesionToadsThinningPass:
    """Run the next ``initThinningLabels``/``propagateThinning`` pair.

    The prior iteration must have continued through standard-mode growth and
    refreshed relations. Its grown segmentation and frontier state are used
    exactly as the selected wrapper's next-loop inputs. Class memberships and
    auxiliary lesion maps are carried from that branch's membership result;
    ``transformed_shape`` must be the already transformed atlas shape maps for
    this iteration's post-registration state. Those maps are fixed only for
    this thinning seam; this adapter does not perform the next registration.
    """
    if not isinstance(previous_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("previous_branch must be a LesionToadsPostRegistrationBranch")
    if previous_branch.stop or previous_branch.relations is None:
        raise ValueError("the prior branch must continue and refresh relations")
    if previous_branch.growing is None:
        raise ValueError("the prior branch must contain its standard-mode growth state")
    membership_result = previous_branch.membership_result
    membership = getattr(membership_result, "membership", None)
    if membership is None:
        raise TypeError("the prior branch must carry its membership result")
    types = tuple(str(value).lower() for value in object_types)
    if len(types) != membership.class_memberships.shape[-1]:
        raise ValueError("object_types must match the membership class count")
    if not 0 <= int(lesion_class) < len(types):
        raise ValueError("lesion_class must index the membership class profile")

    growth = previous_branch.growing
    return propagate_lesion_toads_thinning(
        previous_branch.segmentation,
        membership.class_memberships,
        types,
        membership.lesion_membership,
        membership.black_hole_membership,
        transformed_shape,
        growth.skeleton,
        growth.available,
        growth.ordering,
        lesion_class=int(lesion_class),
        volume_ratio=volume_ratio,
        limit=limit,
        use_lesion_max=use_lesion_max,
    )
