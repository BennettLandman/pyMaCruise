"""Compose the selected LesionTOADS post-registration wrapper branch."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_growth import LesionToadsGrowingPass, propagate_lesion_toads_growing
from .lesion_registration_handoff import LesionRegistrationHandoff
from .lesion_relations import LesionRelationList, lesion_compute_relations
from .lesion_state import LesionToadsFinalClassification, lesion_toads_final_classification


@dataclass(frozen=True)
class LesionToadsPostRegistrationBranch:
    """State and stage trace after one active LesionTOADS alignment update."""

    membership_result: Any
    registration_handoff: LesionRegistrationHandoff
    energy_ratio: np.float32
    convergence_stop: bool
    iteration_limit_stop: bool
    stop: bool
    stop_reason: str | None
    run_growing: bool
    growing: LesionToadsGrowingPass | None
    segmentation: NDArray[np.int16]
    final_classification: LesionToadsFinalClassification | None
    relations: LesionRelationList | None
    stage_order: tuple[str, ...]


def lesion_toads_post_registration_branch(
    membership_result: Any,
    registration_handoff: LesionRegistrationHandoff,
    *,
    previous_energy: float,
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
) -> LesionToadsPostRegistrationBranch:
    """Run only the selected wrapper's branch immediately after alignment.

    Energy and updated memberships come from P39zzzzs, and the already-run
    single-transform alignment state comes from P39zzzzt. The selected wrapper
    recomputes its float32 energy ratio, tests convergence before the iteration
    ceiling, grows only outside ``narrow_band_competition`` when ``!stop`` or
    the original signed ``nIterMax`` is positive, then either transfers final
    labels/subtypes or refreshes relations. The LesionTOADS branch does not
    update base-TOADS boundary-distance or closest-label maps here.
    """
    if not isinstance(registration_handoff, LesionRegistrationHandoff):
        raise TypeError("registration_handoff must be a LesionRegistrationHandoff")
    if not isinstance(algorithm_mode, str):
        raise TypeError("algorithm_mode must be a string")
    if not isinstance(has_lesion, (bool, np.bool_)) or not isinstance(
        use_lesion_max, (bool, np.bool_)
    ):
        raise TypeError("lesion flags must be boolean")
    if isinstance(iteration, (bool, np.bool_)) or not isinstance(
        iteration, (int, np.integer)
    ) or iteration < 1:
        raise ValueError("iteration must be a positive integer")
    if isinstance(n_iter_max, (bool, np.bool_)) or not isinstance(
        n_iter_max, (int, np.integer)
    ) or n_iter_max == 0:
        raise ValueError("n_iter_max must be a nonzero integer")
    iterations = abs(int(n_iter_max))
    if iteration > iterations:
        raise ValueError("iteration cannot exceed abs(n_iter_max)")
    if isinstance(relation_limit, (bool, np.bool_)) or not isinstance(
        relation_limit, (int, np.integer)
    ) or relation_limit < 0:
        raise ValueError("relation_limit must be a nonnegative integer")

    membership = getattr(membership_result, "membership", None)
    if membership is None or not hasattr(membership, "objective"):
        raise TypeError("membership_result must be a P39zzzzs statistics bridge result")
    updated_memberships = np.asarray(membership.class_memberships, dtype=np.float32)
    aligned_memberships = np.asarray(registration_handoff.memberships, dtype=np.float32)
    if updated_memberships.shape != aligned_memberships.shape or not np.array_equal(
        updated_memberships, aligned_memberships
    ):
        raise ValueError("registration handoff must contain the updated P39zzzzs memberships")

    previous = np.float32(previous_energy)
    energy = np.float32(membership.objective)
    distance = np.float32(max_distance)
    if not np.isfinite(previous) or not np.isfinite(energy):
        raise ValueError("previous and current energy must be finite")
    if not np.isfinite(distance) or distance < np.float32(0.0):
        raise ValueError("max_distance must be finite and nonnegative")
    difference = np.float32(previous - energy)
    numerator = np.float32(np.float32(2.0) * difference)
    denominator = np.float32(previous + energy)
    with np.errstate(divide="ignore", invalid="ignore"):
        energy_ratio = np.float32(np.abs(np.float32(numerator / denominator)))
    convergence_stop = bool(energy_ratio < distance)
    iteration_limit_stop = int(iteration) >= iterations
    stop = convergence_stop or iteration_limit_stop
    stop_reason = (
        "convergence_and_iteration_limit"
        if convergence_stop and iteration_limit_stop
        else "convergence"
        if convergence_stop
        else "iteration_limit"
        if iteration_limit_stop
        else None
    )

    narrow_band = algorithm_mode == "narrow_band_competition"
    run_growing = bool((not narrow_band) and ((not stop) or (int(n_iter_max) > 0)))
    stages = [
        "membership_energy",
        "atlas_alignment",
        "energy_ratio",
        "convergence_check",
        "iteration_limit_check",
    ]
    source_segmentation = np.asarray(segmentation)
    if source_segmentation.ndim != 3:
        raise ValueError("segmentation must be a 3-D label volume")
    if source_segmentation.shape != updated_memberships.shape[:3]:
        raise ValueError("segmentation and updated memberships must share a grid")

    growing = None
    if run_growing:
        stages.append("propagateGrowing")
        growing = propagate_lesion_toads_growing(
            source_segmentation,
            skeleton,
            updated_memberships,
            ordering,
            class_volume,
            tuple(str(value) for value in object_types),
            membership.lesion_membership,
            membership.black_hole_membership,
            lesion_class=lesion_class,
            use_lesion_max=bool(use_lesion_max),
        )
        segmentation_after = growing.segmentation.copy()
    else:
        segmentation_after = source_segmentation.astype(np.int16, copy=True)

    final_classification = None
    relations = None
    if stop:
        stages.extend(("updateClassificationFromLabels", "updateLesionClassifcation"))
        final_classification = lesion_toads_final_classification(
            segmentation_after,
            template_labels,
            classification,
            updated_memberships,
            membership.lesion_membership,
            membership.black_hole_membership,
            lesion_class,
            has_lesion=bool(has_lesion),
        )
    else:
        stages.append("computeRelations")
        relations = lesion_compute_relations(
            segmentation_after, updated_memberships.shape[-1], int(relation_limit)
        )

    return LesionToadsPostRegistrationBranch(
        membership_result=membership_result,
        registration_handoff=registration_handoff,
        energy_ratio=energy_ratio,
        convergence_stop=convergence_stop,
        iteration_limit_stop=iteration_limit_stop,
        stop=stop,
        stop_reason=stop_reason,
        run_growing=run_growing,
        growing=growing,
        segmentation=segmentation_after,
        final_classification=final_classification,
        relations=relations,
        stage_order=tuple(stages),
    )
