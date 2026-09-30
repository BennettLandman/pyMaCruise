"""Repeat the selected fixed-input standard LesionTOADS iteration.

This driver starts from a supplied continuing post-registration state and
repeats the already-composed P39v membership update and P39w registration /
stop branch. It does not initialize atlas state or estimate fixed maps.

Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md``, describes
an iterative lesion-aware statistical and topological segmentation. This
driver captures only the selected VandyCRUISE wrapper's repeated transition
and stop order for fixed inputs; these control-flow details come from the
source, not the paper. Synthetic Java-harness comparisons cover the stated
branch/loop cases and do not validate full LesionTOADS convergence or output
equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_continuing_iteration import LesionToadsContinuingIteration, run_lesion_toads_continuing_iteration
from .lesion_continuing_post_membership import run_lesion_toads_continuing_post_membership_branch
from .lesion_post_registration import LesionToadsPostRegistrationBranch


@dataclass(frozen=True)
class LesionToadsLoopPass:
    """One P39v/P39w pass with its one-based selected-wrapper iteration."""

    iteration: int
    continuing_iteration: LesionToadsContinuingIteration
    post_membership_branch: LesionToadsPostRegistrationBranch


@dataclass(frozen=True)
class LesionToadsFixedInputLoop:
    """Repeated selected loop result, ending at a source stop decision."""

    passes: tuple[LesionToadsLoopPass, ...]
    final_branch: LesionToadsPostRegistrationBranch
    stop_reason: str
    stage_order: tuple[str, ...]


def run_lesion_toads_fixed_input_loop(
    initial_branch: LesionToadsPostRegistrationBranch,
    images: ArrayLike,
    transformed_shape: ArrayLike,
    class_shape_values: ArrayLike,
    shape_priors: Sequence[ArrayLike],
    *,
    starting_iteration: int,
    n_iter_max: int,
    max_distance: float,
    relation_limit: int,
    volume_ratio: float,
    thinning_limit: float,
    bridge_options: Mapping[str, Any],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    object_types: Sequence[str],
    classification: ArrayLike,
    template_labels: ArrayLike,
    lesion_class: int,
    has_lesion: bool = True,
    use_lesion_max: bool = False,
    align_shapes: bool = True,
    init_align_iter: int = 10,
    use_lesion_max_thinning: bool = False,
) -> LesionToadsFixedInputLoop:
    """Run P39v then P39w until convergence or the signed iteration ceiling.

    ``initial_branch`` is the selected wrapper state after a completed,
    non-stopping iteration. Images, correction fields, shape values/priors,
    and other maps in ``bridge_options`` stay fixed throughout this bounded
    slice. Memberships, objective, transform, segmentation, growth frontier,
    and relation state advance through each P39v/P39w result.
    """
    if not isinstance(initial_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("initial_branch must be a LesionToadsPostRegistrationBranch")
    if initial_branch.stop or initial_branch.relations is None or initial_branch.growing is None:
        raise ValueError("initial_branch must be continuing with refreshed relations and growth state")
    if isinstance(starting_iteration, (bool, np.bool_)) or not isinstance(
        starting_iteration, (int, np.integer)
    ) or starting_iteration < 1:
        raise ValueError("starting_iteration must be a positive completed iteration count")
    if isinstance(n_iter_max, (bool, np.bool_)) or not isinstance(
        n_iter_max, (int, np.integer)
    ) or n_iter_max == 0:
        raise ValueError("n_iter_max must be a nonzero signed integer")
    iteration_ceiling = abs(int(n_iter_max))
    if int(starting_iteration) >= iteration_ceiling:
        raise ValueError("a continuing initial state must precede the iteration ceiling")
    branch = initial_branch
    loop_passes: list[LesionToadsLoopPass] = []
    stages: list[str] = []
    for iteration in range(int(starting_iteration) + 1, iteration_ceiling + 1):
        continuing = run_lesion_toads_continuing_iteration(
            branch,
            images,
            transformed_shape,
            class_shape_values,
            volume_ratio=volume_ratio,
            thinning_limit=thinning_limit,
            bridge_options=bridge_options,
            use_lesion_max=use_lesion_max_thinning,
        )
        post_branch = run_lesion_toads_continuing_post_membership_branch(
            branch,
            continuing,
            shape_priors,
            previous_energy=float(np.float32(branch.membership_result.membership.objective)),
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            object_types=object_types,
            classification=classification,
            template_labels=template_labels,
            lesion_class=lesion_class,
            iteration=iteration,
            n_iter_max=int(n_iter_max),
            max_distance=max_distance,
            algorithm_mode="standard",
            relation_limit=relation_limit,
            has_lesion=has_lesion,
            use_lesion_max=use_lesion_max,
            align_shapes=align_shapes,
            init_align_iter=init_align_iter,
        )
        loop_passes.append(LesionToadsLoopPass(iteration, continuing, post_branch))
        stages.extend(post_branch.stage_order)
        branch = post_branch
        if branch.stop:
            stop_reason = branch.stop_reason
            break
    else:
        raise RuntimeError("selected wrapper loop reached its ceiling without a stop decision")

    return LesionToadsFixedInputLoop(
        passes=tuple(loop_passes),
        final_branch=branch,
        stop_reason=str(stop_reason),
        stage_order=tuple(stages),
    )
