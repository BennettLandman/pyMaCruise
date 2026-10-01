"""Compose repeated fixed-input LesionTOADS passes through P39v and P39w."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_continuing_iteration import LesionToadsContinuingIteration, run_lesion_toads_continuing_iteration
from .lesion_post_registration import LesionToadsPostRegistrationBranch
from .lesion_registration_branch_composition import (
    LesionToadsRegistrationBranchComposition,
    finish_lesion_toads_iteration_after_membership,
)


@dataclass(frozen=True)
class LesionToadsRepeatedIterationPass:
    """One continuing membership pass and its following registration branch."""

    continuing: LesionToadsContinuingIteration
    registration_branch: LesionToadsRegistrationBranchComposition
    previous_energy: np.float32
    initial_transform: np.ndarray
    starting_class_memberships: np.ndarray
    starting_segmentation: np.ndarray
    starting_relation_map: np.ndarray
    starting_lesion_membership: np.ndarray
    starting_black_hole_membership: np.ndarray
    starting_skeleton: np.ndarray
    starting_available: np.ndarray
    starting_ordering: np.ndarray
    starting_class_volume: np.ndarray


@dataclass(frozen=True)
class LesionToadsRepeatedIterationResult:
    """Selected fixed-input loop result through its first stop condition."""

    passes: tuple[LesionToadsRepeatedIterationPass, ...]
    final_branch: LesionToadsPostRegistrationBranch
    stop_reason: str
    stage_order: tuple[tuple[str, ...], ...]


def run_lesion_toads_repeated_iterations(
    initial_branch: LesionToadsPostRegistrationBranch,
    images: ArrayLike,
    transformed_shape: ArrayLike,
    class_shape_values: ArrayLike,
    *,
    shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    bridge_options: Mapping[str, Any],
    volume_ratio: float,
    thinning_limit: float,
    classification: ArrayLike,
    template_labels: ArrayLike,
    object_types: Sequence[str],
    lesion_class: int,
    n_iter_max: int,
    max_distance: float,
    algorithm_mode: str,
    relation_limit: int,
    align_shapes: bool = True,
    init_align_iter: int = 10,
    use_lesion_max: bool = False,
) -> LesionToadsRepeatedIterationResult:
    """Run selected continuing iterations until convergence or signed limit.

    Geometry, atlas/shape maps, correction fields, and other `bridge_options`
    are fixed across passes. Membership, topology, relation, objective, and
    transform states flow from each P39v/P39w result into the next pass. This
    is a bounded standard input-state loop; initialization, field estimation,
    atlas resampling, and external workflow outputs remain outside it.
    """
    if not isinstance(initial_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("initial_branch must be a LesionToadsPostRegistrationBranch")
    if initial_branch.stop or initial_branch.relations is None or initial_branch.growing is None:
        raise ValueError("initial branch must be a continuing post-registration state")
    if isinstance(n_iter_max, (bool, np.bool_)) or not isinstance(n_iter_max, (int, np.integer)) or n_iter_max == 0:
        raise ValueError("n_iter_max must be a nonzero integer")
    iterations = abs(int(n_iter_max))

    branch = initial_branch
    transform = np.asarray(initial_transform, dtype=np.float32).copy()
    passes: list[LesionToadsRepeatedIterationPass] = []
    for iteration in range(1, iterations + 1):
        prior_membership = branch.membership_result.membership
        assert branch.relations is not None
        starting_class_memberships = prior_membership.class_memberships.copy()
        starting_segmentation = branch.segmentation.copy()
        starting_relation_map = branch.relations.relation_map.copy()
        starting_lesion_membership = prior_membership.lesion_membership.copy()
        starting_black_hole_membership = prior_membership.black_hole_membership.copy()
        assert branch.growing is not None
        starting_skeleton = branch.growing.skeleton.copy()
        starting_available = branch.growing.available.copy()
        starting_ordering = branch.growing.ordering.copy()
        starting_class_volume = branch.growing.class_volume.copy()
        continuing = run_lesion_toads_continuing_iteration(
            branch,
            images,
            transformed_shape,
            class_shape_values,
            volume_ratio=volume_ratio,
            thinning_limit=thinning_limit,
            bridge_options=bridge_options,
            use_lesion_max=use_lesion_max,
        )
        previous_energy = np.float32(branch.membership_result.membership.objective)
        transform_before_registration = transform.copy()
        registration_branch = finish_lesion_toads_iteration_after_membership(
            continuing,
            previous_energy=float(previous_energy),
            shape_priors=shape_priors,
            initial_transform=transform,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            classification=classification,
            template_labels=template_labels,
            object_types=object_types,
            lesion_class=lesion_class,
            iteration=iteration,
            n_iter_max=int(n_iter_max),
            max_distance=max_distance,
            algorithm_mode=algorithm_mode,
            relation_limit=relation_limit,
            align_shapes=align_shapes,
            init_align_iter=init_align_iter,
        )
        passes.append(
            LesionToadsRepeatedIterationPass(
                continuing,
                registration_branch,
                previous_energy,
                transform_before_registration,
                starting_class_memberships,
                starting_segmentation,
                starting_relation_map,
                starting_lesion_membership,
                starting_black_hole_membership,
                starting_skeleton,
                starting_available,
                starting_ordering,
                starting_class_volume,
            )
        )
        branch = registration_branch.branch
        transform = registration_branch.registration.transform.copy()
        if branch.stop:
            assert branch.stop_reason is not None
            return LesionToadsRepeatedIterationResult(
                tuple(passes), branch, branch.stop_reason,
                tuple(value.registration_branch.stage_order for value in passes),
            )

    raise RuntimeError("selected iteration ceiling completed without a stopping branch")
