"""Drive the selected 13-class LesionTOADS registered loop.

This joins already source-checked transition, registration, and wrapper
branch adapters. Atlas shapes and, when atlas names are provided, distance
maps are refreshed after each thinning pass. Correction fields remain caller
inputs; preprocessing, initial registration, and complete workflow I/O are
outside this runner.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence
from time import perf_counter

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
from .lesion_registered_initial_membership import LesionInitialDistanceMaps
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
    correction_polynomial_degree: int | None = None
    initial_correction_degree: int = 0
    align_shapes: bool = True
    init_align_iter: int = 10
    has_lesion: bool = True
    use_lesion_max: bool = False
    retain_growth_event_states: bool = True
    retain_iteration_details: bool = False
    atlas_names: Sequence[str] | None = None
    maximum_ventricle_distance: int = 2
    maximum_gm_distance: int = 3
    maximum_brainstem_distance: int = 25
    initial_distance_maps: LesionInitialDistanceMaps | None = None


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
    refreshed_distance_maps: LesionInitialDistanceMaps | None = None
    correction_degree: int = 0


@dataclass(frozen=True)
class LesionToadsDescriptorLoopResult:
    """Registered iterations ending at a source stop or requested checkpoint."""

    passes: tuple[LesionToadsDescriptorLoopPass, ...]
    final_branch: LesionToadsPostRegistrationBranch
    stop_reason: str
    stage_order: tuple[tuple[str, ...], ...]
    final_class_centroids: np.ndarray
    final_lesion_centroids: np.ndarray
    final_black_hole_centroids: np.ndarray
    final_distance_maps: LesionInitialDistanceMaps | None = None
    final_fields: np.ndarray | None = None
    final_correction_degree: int = 0


def run_lesion_toads_descriptor_loop(
    initial_branch: LesionToadsPostRegistrationBranch,
    initial_registration,
    *,
    starting_iteration: int,
    options: LesionToadsDescriptorLoopOptions,
    max_passes: int | None = None,
    progress_callback: Callable[[dict[str, Any]], None] | None = None,
    stage_progress_callback: Callable[[dict[str, Any]], None] | None = None,
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
    if progress_callback is not None and not callable(progress_callback):
        raise TypeError("progress_callback must be callable when supplied")
    if stage_progress_callback is not None and not callable(stage_progress_callback):
        raise TypeError("stage_progress_callback must be callable when supplied")
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
    if max_passes is not None and (
        isinstance(max_passes, (bool, np.bool_))
        or not isinstance(max_passes, (int, np.integer))
        or max_passes < 1
    ):
        raise ValueError("max_passes must be a positive integer when supplied")
    if not isinstance(options.source_zero_variance_nan, (bool, np.bool_)):
        raise TypeError("source_zero_variance_nan must be boolean")
    if not isinstance(options.retain_iteration_details, (bool, np.bool_)):
        raise TypeError("retain_iteration_details must be boolean")
    if not isinstance(options.retain_growth_event_states, (bool, np.bool_)):
        raise TypeError("retain_growth_event_states must be boolean")
    if options.correction_polynomial_degree is not None and (
        isinstance(options.correction_polynomial_degree, (bool, np.bool_))
        or not isinstance(options.correction_polynomial_degree, (int, np.integer))
        or int(options.correction_polynomial_degree) not in (1, 2, 3)
    ):
        raise ValueError("correction_polynomial_degree must be one of 1, 2, 3, or None")
    if (
        isinstance(options.initial_correction_degree, (bool, np.bool_))
        or not isinstance(options.initial_correction_degree, (int, np.integer))
        or options.initial_correction_degree < 0
        or (
            options.correction_polynomial_degree is not None
            and options.initial_correction_degree > options.correction_polynomial_degree
        )
    ):
        raise ValueError("initial_correction_degree must be nonnegative and no greater than the fit ceiling")
    images = np.asarray(options.images, dtype=np.float32)
    if images.ndim not in (3, 4):
        raise ValueError("images must be a 3-D volume or channel-last 4-D array")
    if options.fields is not None:
        supplied_fields = np.asarray(options.fields, dtype=np.float32)
        if supplied_fields.shape != images.shape or not np.isfinite(supplied_fields).all():
            raise ValueError("fields must be finite and match the loop image grid and channels")
    elif options.correction_polynomial_degree is not None:
        raise ValueError("fields are required when correction-field refitting is enabled")

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
    ventricle_distances = np.asarray(options.ventricle_distances)
    gray_matter_distances = np.asarray(options.gray_matter_distances)
    brainstem_distances = np.asarray(options.brainstem_distances)
    final_distance_maps = None
    current_fields = (
        None
        if options.fields is None
        else np.array(options.fields, dtype=np.float32, copy=True)
    )
    correction_degree = int(options.initial_correction_degree)
    passes: list[LesionToadsDescriptorLoopPass] = []
    loop_started = perf_counter()

    last_iteration = maximum if max_passes is None else min(
        maximum, int(starting_iteration) + int(max_passes)
    )
    for iteration in range(int(starting_iteration) + 1, last_iteration + 1):
        phase_started = perf_counter()
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "registered_shape_priors",
                "scope": "pass",
                "status": "started",
            })
        shape_values = transform_registered_lesion_shape_priors(
            shape_registration,
            options.atlas_shape_priors,
            image_dimensions_xyz=tuple(int(value) for value in images.shape[:3]),
            image_resolutions_mm=options.image_resolutions_mm,
            template_resolutions_mm=options.template_resolutions_mm,
        )
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "registered_shape_priors",
                "scope": "pass",
                "status": "completed",
                "elapsed_seconds": float(perf_counter() - phase_started),
            })
        relation_state = branch.relations
        assert relation_state is not None
        relation_table = relation_state.relation_membership[relation_state.relation_map]
        phase_started = perf_counter()
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "continuing_segmentation_statistics_memberships",
                "scope": "pass",
                "status": "started",
            })
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
            ventricle_distances=ventricle_distances,
            gray_matter_distances=gray_matter_distances,
            brainstem_distances=brainstem_distances,
            centroid_options=options.centroid_options,
            channel_weight_options=options.channel_weight_options,
            outlier_energy=options.outlier_energy,
            lesion_power_exponent=options.lesion_power_exponent,
            fields=current_fields,
            flair_channel=options.flair_channel,
            spread=options.spread,
            source_zero_variance_nan=options.source_zero_variance_nan,
            use_lesion_max=options.use_lesion_max,
            atlas_names=options.atlas_names,
            maximum_ventricle_distance=options.maximum_ventricle_distance,
            maximum_gm_distance=options.maximum_gm_distance,
            maximum_brainstem_distance=options.maximum_brainstem_distance,
            has_lesion=options.has_lesion,
            previous_distance_maps=final_distance_maps or options.initial_distance_maps,
            correction_polynomial_degree=options.correction_polynomial_degree,
            correction_degree=correction_degree,
            correction_max_distance=options.max_distance,
            stage_progress_callback=(
                None
                if stage_progress_callback is None
                else lambda event: stage_progress_callback({
                    "iteration": int(iteration),
                    "scope": "continuing",
                    **event,
                })
            ),
        )
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "continuing_segmentation_statistics_memberships",
                "scope": "pass",
                "status": "completed",
                "elapsed_seconds": float(perf_counter() - phase_started),
            })
        if transition.inhomogeneity is not None:
            correction_degree = transition.inhomogeneity.degree
            current_fields = np.stack(transition.inhomogeneity.fields, axis=-1)
        if transition.refreshed_distance_maps is not None:
            final_distance_maps = transition.refreshed_distance_maps
            ventricle_distances = final_distance_maps.ventricle
            gray_matter_distances = final_distance_maps.gray_matter
            brainstem_distances = final_distance_maps.brainstem
        elif final_distance_maps is None:
            final_distance_maps = options.initial_distance_maps
        phase_started = perf_counter()
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "registration_and_branching",
                "scope": "pass",
                "status": "started",
            })
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
            retain_growth_event_states=options.retain_growth_event_states,
        )
        if stage_progress_callback is not None:
            stage_progress_callback({
                "iteration": int(iteration),
                "stage": "registration_and_branching",
                "scope": "pass",
                "status": "completed",
                "elapsed_seconds": float(perf_counter() - phase_started),
            })
        passes.append(
            LesionToadsDescriptorLoopPass(
                iteration=iteration,
                energy_ratio=np.float32(transition.membership.energy_ratio.energy_ratio),
                registration_transform=registered.registration.registration.transform.copy(),
                stop_reason=branch.stop_reason,
                stage_order=(*transition.stage_order, *registered.stage_order),
                transformed_shape_values=(
                    shape_values.copy() if options.retain_iteration_details else None
                ),
                transition=transition if options.retain_iteration_details else None,
                registered=registered if options.retain_iteration_details else None,
                refreshed_distance_maps=transition.refreshed_distance_maps,
                correction_degree=correction_degree,
            )
        )
        branch = registered.control.branch
        transform = registered.registration.registration.transform.copy()
        shape_registration = registered.registration.registration
        class_centroids = transition.statistics.statistics.centroids.class_centroids
        lesion_centroids = transition.statistics.statistics.centroids.lesion_centroids
        black_hole_centroids = transition.statistics.statistics.centroids.black_hole_centroids
        if progress_callback is not None:
            progress_callback({
                "iteration": int(iteration),
                "energy_ratio": float(transition.membership.energy_ratio.energy_ratio),
                "stop": bool(branch.stop),
                "stop_reason": branch.stop_reason,
                "stage_order": list(passes[-1].stage_order),
                "elapsed_seconds": float(perf_counter() - loop_started),
            })
        if branch.stop:
            assert branch.stop_reason is not None
            return LesionToadsDescriptorLoopResult(
                tuple(passes),
                branch,
                branch.stop_reason,
                tuple(item.stage_order for item in passes),
                class_centroids.copy(),
                lesion_centroids.copy(),
                black_hole_centroids.copy(),
                final_distance_maps,
                None if current_fields is None else current_fields.copy(),
                correction_degree,
            )

    if max_passes is not None:
        return LesionToadsDescriptorLoopResult(
            tuple(passes),
            branch,
            "checkpoint_boundary",
            tuple(item.stage_order for item in passes),
            class_centroids.copy(),
            lesion_centroids.copy(),
            black_hole_centroids.copy(),
            final_distance_maps,
            None if current_fields is None else current_fields.copy(),
            correction_degree,
        )

    raise RuntimeError("selected iteration ceiling completed without a stopping branch")
