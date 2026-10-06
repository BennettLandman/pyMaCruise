"""Continue a selected LesionTOADS update after atlas registration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from numpy.typing import ArrayLike

from .lesion_next_iteration_registration import LesionNextIterationRegistration
from .lesion_post_registration import (
    LesionToadsPostRegistrationBranch,
    lesion_toads_post_registration_branch,
)


@dataclass(frozen=True)
class LesionNextIterationBranch:
    """Registered update plus the wrapper's convergence/growth/output branch."""

    registration: LesionNextIterationRegistration
    branch: LesionToadsPostRegistrationBranch
    stage_order: tuple[str, ...]


def continue_after_next_iteration_registration(
    registration: LesionNextIterationRegistration,
    *,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volume: ArrayLike,
    object_types: Sequence[str],
    classification: ArrayLike,
    template_labels: ArrayLike,
    lesion_class: int,
    iteration: int,
    n_iter_max: int,
    max_distance: float,
    algorithm_mode: str,
    relation_limit: int,
    has_lesion: bool = True,
    use_lesion_max: bool = False,
    progress_callback: Callable[[str], None] | None = None,
) -> LesionNextIterationBranch:
    """Apply selected post-registration stop, optional growth, and refresh.

    The membership pass's source-computed ratio is carried into the wrapper
    branch. Registration has already completed, matching the selected order in
    ``AlgorithmLesionToads`` where ``registerShapes`` precedes convergence.
    This composes one supplied state only; it does not start another iteration.
    """
    if not isinstance(registration, LesionNextIterationRegistration):
        raise TypeError("registration must be a LesionNextIterationRegistration")
    update = registration.membership_update
    branch = lesion_toads_post_registration_branch(
        update,
        registration.registration,
        previous_energy=float(update.energy_ratio.previous_energy),
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volume=class_volume,
        object_types=object_types,
        classification=classification,
        template_labels=template_labels,
        lesion_class=lesion_class,
        iteration=iteration,
        n_iter_max=n_iter_max,
        max_distance=max_distance,
        algorithm_mode=algorithm_mode,
        relation_limit=relation_limit,
        has_lesion=has_lesion,
        use_lesion_max=use_lesion_max,
        source_energy_ratio=update.energy_ratio,
        progress_callback=progress_callback,
    )
    return LesionNextIterationBranch(
        registration,
        branch,
        registration.stage_order + branch.stage_order[3:],
    )
