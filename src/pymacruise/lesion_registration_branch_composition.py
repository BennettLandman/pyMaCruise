"""Connect one P39v result to the selected LesionTOADS alignment branch.

This composes existing selected-source kernels; it does not add registration,
membership, growth, or relation mathematics of its own.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from numpy.typing import ArrayLike

from .lesion_continuing_iteration import LesionToadsContinuingIteration
from .lesion_post_registration import LesionToadsPostRegistrationBranch, lesion_toads_post_registration_branch
from .lesion_registration_handoff import LesionRegistrationHandoff, register_lesion_membership_handoff


@dataclass(frozen=True)
class LesionToadsRegistrationBranchComposition:
    """P39v membership output followed by registration and wrapper branching."""

    registration: LesionRegistrationHandoff
    branch: LesionToadsPostRegistrationBranch
    stage_order: tuple[str, ...]


def finish_lesion_toads_iteration_after_membership(
    continuing: LesionToadsContinuingIteration,
    *,
    previous_energy: float,
    shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    classification: ArrayLike,
    template_labels: ArrayLike,
    object_types: Sequence[str],
    lesion_class: int,
    iteration: int,
    n_iter_max: int,
    max_distance: float,
    algorithm_mode: str,
    relation_limit: int,
    align_shapes: bool = True,
    init_align_iter: int = 10,
    has_lesion: bool = True,
    use_lesion_max: bool = False,
) -> LesionToadsRegistrationBranchComposition:
    """Run the selected registration then convergence/growth/output branch.

    The supplied result is P39v's just-computed membership/objective state.
    Shape/spatial maps and geometry are fixed call inputs. The previous
    iteration state supplies the topology/growth scaffolding captured by P39v
    thinning; no further thinning or membership pass is run here.
    """
    if not isinstance(continuing, LesionToadsContinuingIteration):
        raise TypeError("continuing must be a LesionToadsContinuingIteration")
    membership_result = continuing.statistics_energy
    thinning = continuing.thinning.run
    registration = register_lesion_membership_handoff(
        membership_result,
        shape_priors,
        initial_transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    branch = lesion_toads_post_registration_branch(
        membership_result,
        registration,
        previous_energy=previous_energy,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=thinning.ordering,
        class_volume=thinning.class_volume,
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
    )
    # In AlgorithmLesionToads.runAlgorithm(), the ratio is set directly after
    # computeMemberships and registration follows it. The post-registration
    # helper recalculates the ratio only to expose the branch predicate.
    # Collapse P39v's separately tested helper kernels to the selected Java
    # wrapper calls; computeMemberships() encloses its internal substeps.
    stage_order = (
        "propagateThinning",
        "computeCentroids",
        "computeVariances",
        "computeChannelWeights",
        "computeMemberships",
        "energyRatio",
        "registerShapes",
        "convergenceCheck",
        "iterationLimitCheck",
        *(("propagateGrowing",) if branch.run_growing else ()),
        *(("updateClassificationFromLabels", "updateLesionClassifcation") if branch.stop else ("computeRelations",)),
    )
    return LesionToadsRegistrationBranchComposition(registration, branch, stage_order)
