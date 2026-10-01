"""Compose initial LesionTOADS memberships into the first thinning pass.

The selected wrapper computes initial memberships with ``GMflag=false``,
then sets the flag and enters the main loop. This adapter carries those
initial membership maps directly into ``initThinningLabels`` and
``propagateThinning``. It deliberately stops after this first pass; the later
``GMflag=true`` membership update remains a separate iteration stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from numpy.typing import ArrayLike

from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_initial_loop_entry import LesionToadsInitialLoopEntry


@dataclass(frozen=True)
class LesionInitialThinningComposition:
    """Selected loop-entry state and its first source-ordered thinning pass."""

    loop_entry: LesionToadsInitialLoopEntry
    thinning: LesionToadsThinningPass
    stage_order: tuple[str, ...]


def run_lesion_initial_thinning_pass(
    loop_entry: LesionToadsInitialLoopEntry,
    segmentation: ArrayLike,
    object_types: Sequence[str],
    transformed_shape: ArrayLike,
    initial_skeleton: ArrayLike,
    initial_available: ArrayLike,
    initial_ordering: ArrayLike,
    *,
    lesion_class: int,
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
    empty_label: int = -1,
) -> LesionInitialThinningComposition:
    """Run first-loop thinning on the exact initial membership result.

    ``segmentation``, transformed atlas shapes, and the initial topology/frontier
    state are fixed caller inputs. No later membership computation is run.
    """
    if not isinstance(loop_entry, LesionToadsInitialLoopEntry):
        raise TypeError("loop_entry must be a LesionToadsInitialLoopEntry")
    if not loop_entry.enter_loop or loop_entry.next_stage != "propagateThinning":
        raise ValueError("selected loop entry does not reach propagateThinning")
    if not loop_entry.gm_flag:
        raise ValueError("selected wrapper sets GMflag before entering thinning")
    initial = loop_entry.initial_membership
    thinning = propagate_lesion_toads_thinning(
        segmentation,
        initial.class_memberships,
        object_types,
        initial.lesion_membership,
        initial.black_hole_membership,
        transformed_shape,
        initial_skeleton,
        initial_available,
        initial_ordering,
        lesion_class=lesion_class,
        volume_ratio=volume_ratio,
        limit=limit,
        use_lesion_max=use_lesion_max,
        empty_label=empty_label,
    )
    return LesionInitialThinningComposition(
        loop_entry,
        thinning,
        loop_entry.stage_order + ("propagateThinning:initThinningLabels", "propagateThinning:heapLoop"),
    )
