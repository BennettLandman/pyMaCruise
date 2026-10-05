"""Drive the selected 13-class LesionTOADS registered loop.

This joins already source-checked transition, registration, and wrapper
branch adapters. Atlas shapes are resampled after each registration. The
selected distance maps and correction fields are supplied by the caller and
remain fixed across this bounded loop; preprocessing, initial registration,
and complete workflow I/O are outside this runner.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_descriptor_continuing_iteration import (
    LesionToadsDescriptorContinuingIteration,
    run_lesion_toads_descriptor_continuing_iteration,
)
from .lesion_descriptor_registered_iteration import (
    LesionToadsDescriptorRegisteredIteration,
    register_and_branch_lesion_toads_descriptor_iteration,
)
from .lesion_post_registration import LesionToadsPostRegistrationBranch
from .lesion_registered_shapes import transform_registered_lesion_shape_priors


@dataclass(frozen=True)
class LesionToadsDescriptorLoopOptions:
    """Fixed inputs and selected wrapper settings for registered iterations."""

    images: ArrayLike
    atlas_shape_priors: Sequence[ArrayLike]
    shape_priors: ArrayLike
    smoothing_priors: ArrayLike
    object_types: Sequence[str]
    lesion_class: int
    volume_ratio: float
    thinning_limit: float
    image_resolutions_mm: tuple[float, float, float]
    template_resolutions_mm: tuple[float, float, float]
    classification: ArrayLike
    template_labels: ArrayLike
    n_iter_max: int
    max_distance: float
    algorithm_mode: str
    relation_limit: int
    initial_class_centroids: ArrayLike
    initial_lesion_centroids: ArrayLike
    initial_black_hole_centroids: ArrayLike
    ventricle_distances: ArrayLike
    gray_matter_distances: ArrayLike
    brainstem_distances: ArrayLike
    centroid_options: Mapping[str, Any] | None
    channel_weight_options: Mapping[str, Any]
    outlier_energy: float
    lesion_power_exponent: float
    fields: ArrayLike | None = None
    flair_channel: int | None = None
    spread: float = 0.5
    source_zero_variance_nan: bool = False
    align_shapes: bool = True
    init_align_iter: int = 10
    has_lesion: bool = True
    use_lesion_max: bool = False
    retain_iteration_details: bool = False


@dataclass(frozen=True)
class LesionToadsDescriptorLoopPass:
    """Compact pass summary, with full arrays retained only on request."""

    iteration: int
    energy_ratio: np.float32
    registration_transform: np.ndarray
    stop_reason: str | None
    stage_order: tuple[str, ...]
    transformed_shape_values: np.ndarray | None = None
    transition: LesionToadsDescriptorContinuingIteration | None = None
    registered: LesionToadsDescriptorRegisteredIteration | None = None


@dataclass(frozen=True)
class LesionToadsDescriptorLoopResult:
    """Registered iterations through the first selected stop condition."""

    passes: tuple[LesionToadsDescriptorLoopPass, ...]
    final_branch: LesionToadsPostRegistrationBranch
    stop_reason: str
    stage_order: tuple[tuple[str, ...], ...]


def run_lesion_toads_descriptor_loop(
    initial_branch: LesionToadsPostRegistrationBranch,
    initial_registration,
    *,
    starting_iteration: int,
    options: LesionToadsDescriptorLoopOptions,
) -> LesionToadsDescriptorLoopResult:
    """Run registered 13-class updates until source convergence or loop cap.

    ``initial_branch`` must be a continuing branch immediately before the
    next thinning/statistics/membership pass. ``starting_iteration`` is the
    last completed source iteration, so the first pass here uses its successor.
    Atlas shape values are recalculated from template priors at each current
    registration transform. Relations, memberships, segmentation, growth
    frontier, and statistics flow from the preceding branch/update.
    """
    if not isinstance(initial_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("initial_branch must be a LesionToadsPostRegistrationBranch")
    if initial_branch.stop or initial_branch.relations is None or initial_branch.growing is None:
        raise ValueError("initial_branch must be a continuing registered state")
    if not isinstance(options, LesionToadsDescriptorLoopOptions):
        raise TypeError("options must be LesionToadsDescriptorLoopOptions")
    if isinstance(starting_iteration, (bool, np.bool_)) or not isinstance(
        starting_iteration, (int, np.integer)
    ) or starting_iteration < 1:
        raise ValueError("starting_iteration must be a positive completed iteration")
    if isinstance(options.n_iter_max, (bool, np.bool_)) or not isinstance(
        options.n_iter_max, (int, np.integer)
    ) or options.n_iter_max == 0:
        raise ValueError("n_iter_max must be a nonzero integer")
    maximum = abs(int(options.n_iter_max))
    if starting_iteration >= maximum:
        raise ValueError("starting_iteration must be below abs(n_iter_max)")
    if not isinstance(options.source_zero_variance_nan, (bool, np.bool_)):
        raise TypeError("source_zero_variance_nan must be boolean")
    if not isinstance(options.retain_iteration_details, (bool, np.bool_)):
        raise TypeError("retain_iteration_details must be boolean")

    images = np.asarray(options.images, dtype=np.float32)
    if images.ndim not in (3, 4):
        raise ValueError("images must be a 3-D volume or channel-last 4-D array")
    transform = np.asarray(getattr(initial_registration, "transform", None), dtype=np.float32).copy()
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_registration must provide a finite six-value rigid transform")
    if getattr(initial_registration, "shape_bounds", None) is None:
        raise ValueError("initial_registration must provide transformed atlas shape bounds")
    shape_registration = initial_registration

    class_centroids = np.asarray(options.initial_class_centroids, dtype=np.float32)
    lesion_centroids = np.asarray(options.initial_lesion_centroids, dtype=np.float32)
    black_hole_centroids = np.asarray(options.initial_black_hole_centroids, dtype=np.float32)
    branch = initial_branch
    passes: list[LesionToadsDescriptorLoopPass] = []

    for iteration in range(int(starting_iteration) + 1, maximum + 1):
        shape_values = transform_registered_lesion_shape_priors(
            shape_registration,
            options.atlas_shape_priors,
            image_dimensions_xyz=tuple(int(value) for value in images.shape[:3]),
            image_resolutions_mm=options.image_resolutions_mm,
            template_resolutions_mm=options.template_resolutions_mm,
        )
        relation_state = branch.relations
        assert relation_state is not None
        relation_table = relation_state.relation_membership[relation_state.relation_map]
        transition = run_lesion_toads_descriptor_continuing_iteration(
            branch,
            images,
            shape_values,
            object_types=options.object_types,
            lesion_class=options.lesion_class,
            volume_ratio=options.volume_ratio,
            thinning_limit=options.thinning_limit,
            initial_class_centroids=class_centroids,
            initial_lesion_centroids=lesion_centroids,
            initial_black_hole_centroids=black_hole_centroids,
            shape_priors=options.shape_priors,
            smoothing_priors=options.smoothing_priors,
            direct_relations=relation_table,
            ventricle_distances=options.ventricle_distances,
            gray_matter_distances=options.gray_matter_distances,
            brainstem_distances=options.brainstem_distances,
            centroid_options=options.centroid_options,
            channel_weight_options=options.channel_weight_options,
            outlier_energy=options.outlier_energy,
            lesion_power_exponent=options.lesion_power_exponent,
            fields=options.fields,
            flair_channel=options.flair_channel,
            spread=options.spread,
            source_zero_variance_nan=options.source_zero_variance_nan,
            use_lesion_max=options.use_lesion_max,
        )
        registered = register_and_branch_lesion_toads_descriptor_iteration(
            transition,
            options.atlas_shape_priors,
            transform,
            image_resolutions_mm=options.image_resolutions_mm,
            template_resolutions_mm=options.template_resolutions_mm,
            object_types=options.object_types,
            classification=options.classification,
            template_labels=options.template_labels,
            lesion_class=options.lesion_class,
            iteration=iteration,
            n_iter_max=int(options.n_iter_max),
            max_distance=options.max_distance,
            algorithm_mode=options.algorithm_mode,
            relation_limit=options.relation_limit,
            align_shapes=options.align_shapes,
            init_align_iter=options.init_align_iter,
            has_lesion=options.has_lesion,
            use_lesion_max=options.use_lesion_max,
        )
        passes.append(
            LesionToadsDescriptorLoopPass(
                iteration=iteration,
                energy_ratio=np.float32(transition.membership.energy_ratio.energy_ratio),
                registration_transform=registered.registration.registration.transform.copy(),
                stop_reason=branch.stop_reason,
                stage_order=registered.stage_order,
                transformed_shape_values=(
                    shape_values.copy() if options.retain_iteration_details else None
                ),
                transition=transition if options.retain_iteration_details else None,
                registered=registered if options.retain_iteration_details else None,
            )
        )
        branch = registered.control.branch
        transform = registered.registration.registration.transform.copy()
        shape_registration = registered.registration.registration
        class_centroids = transition.statistics.statistics.centroids.class_centroids
        lesion_centroids = transition.statistics.statistics.centroids.lesion_centroids
        black_hole_centroids = transition.statistics.statistics.centroids.black_hole_centroids
        if branch.stop:
            assert branch.stop_reason is not None
            return LesionToadsDescriptorLoopResult(
                tuple(passes),
                branch,
                branch.stop_reason,
                tuple(item.stage_order for item in passes),
            )

    raise RuntimeError("selected iteration ceiling completed without a stopping branch")
