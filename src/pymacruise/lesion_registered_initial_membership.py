"""Compose the first LesionTOADS membership pass after registered recropping.

The selected wrapper initializes class centroids before atlas registration,
retains them across the preprocessing recrop, creates zero membership maps in
the segmentation constructor, computes relations/statistics/channel weights,
and runs the first membership update with GMflag already set. This adapter
joins those source-checked operations on a small synthetic or caller-supplied
RAS grid; it does not estimate atlas registration or correction fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_brainstem import lesion_toads_brainstem_distance
from .lesion_gm_distance import lesion_toads_gm_distance
from .lesion_initial_loop_entry import (
    LesionToadsInitialLoopEntry,
    seed_lesion_toads_initial_loop,
)
from .lesion_initial_thinning_composition import (
    LesionInitialThinningComposition,
    run_lesion_initial_thinning_pass,
)
from .lesion_relations import LesionRelationList, lesion_compute_relations
from .lesion_statistics import (
    LesionToadsVariances,
    compute_lesion_toads_channel_weights,
    compute_lesion_toads_variances,
)
from .lesion_standard13_voxel import (
    STANDARD13_NAMES,
    STANDARD13_TYPES,
)
from .lesion_standard13_volume import (
    Standard13VolumeResult,
    run_lesion_toads_standard13_volume,
)
from .lesion_ventricle import lesion_toads_ventricle_distance


@dataclass(frozen=True)
class LesionInitialDistanceMaps:
    ventricle_mask: NDArray[np.bool_]
    ventricle: NDArray[np.int8]
    gray_matter_mask: NDArray[np.bool_]
    gray_matter: NDArray[np.int8]
    brainstem_mask: NDArray[np.bool_]
    brainstem: NDArray[np.int8]


@dataclass(frozen=True)
class LesionRegisteredInitialMembership:
    """First GM-enabled membership pass and the loop's first thinning handoff."""

    relations: LesionRelationList
    variances: LesionToadsVariances
    channel_weights: NDArray[np.float32]
    distance_maps: LesionInitialDistanceMaps
    initial_class_memberships: NDArray[np.float32]
    initial_lesion_membership: NDArray[np.float32]
    initial_black_hole_membership: NDArray[np.float32]
    initial_outlier_membership: NDArray[np.float32]
    initial_lesion_relation: NDArray[np.float32]
    membership: Standard13VolumeResult
    loop_entry: LesionToadsInitialLoopEntry
    first_thinning: LesionInitialThinningComposition
    stage_order: tuple[str, ...]
    excluded_wrapper_stages: tuple[str, ...]


def run_lesion_toads_registered_initial_membership(
    images: ArrayLike,
    segmentation: ArrayLike,
    class_centroids: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    transformed_shapes: ArrayLike,
    skeleton: ArrayLike,
    available: ArrayLike,
    ordering: ArrayLike,
    *,
    relation_limit: int = 3,
    maximum_gm_distance: int = 3,
    maximum_ventricle_distance: int = 2,
    maximum_brainstem_distance: int = 25,
    outlier_energy: float = 0.1,
    lesion_power_exponent: float = 1.2,
    spread: float = 1.0,
    norm_type: str = "Identity",
    iterations: int = 99,
    first_thinning_limit: float = 3.0,
    first_thinning_volume_ratio: float = 0.5,
    source_zero_variance_nan: bool = True,
) -> LesionRegisteredInitialMembership:
    """Run relations, initial statistics, GM-enabled memberships, then thin.

    ``images`` is channel-last ``(x, y, z, 2)`` in T1_MPRAGE/FLAIR order;
    all arrays are on the same already-RAS recropped grid. Class centroids are
    retained pre-registration values. The constructor's current membership
    maps are initialized here to the selected Java defaults (class/lesion/BH
    zero, outlier 0.5, lesion relation 1). Transformed shape maps are supplied
    from the preceding registered-recrop stage.

    The source wrapper calls ``setGMflag`` before its initial
    ``computeMemberships`` (and again afterward), so this first volume pass
    uses ``gm_flag=True``. ``source_zero_variance_nan`` defaults to the
    selected Java IEEE behavior when a computed standard deviation is zero;
    false is available only for a deliberately finite-stabilized comparison.
    """
    data = np.asarray(images, dtype=np.float32)
    labels = np.asarray(segmentation)
    centers = np.asarray(class_centroids, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    black_center = np.asarray(black_hole_centroid, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    skeleton_values = np.asarray(skeleton)
    available_values = np.asarray(available, dtype=np.bool_)
    ordering_values = np.asarray(ordering, dtype=np.float32)
    if data.ndim != 4 or data.shape[-1] != 2 or min(data.shape[:3]) < 3:
        raise ValueError("images must have shape (x>=3,y>=3,z>=3,2) in T1_MPRAGE/FLAIR order")
    spatial = data.shape[:3]
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be integer labels on the image grid")
    if shapes.shape != (*spatial, len(STANDARD13_NAMES)):
        raise ValueError("transformed_shapes must match the image grid and 13-class profile")
    if centers.shape != (2, len(STANDARD13_NAMES)):
        raise ValueError("class_centroids must have shape (2, 13)")
    if lesion_center.shape != (2,) or black_center.shape != (2,):
        raise ValueError("lesion and black-hole centroids must each contain T1 and FLAIR values")
    if skeleton_values.shape != spatial or not np.issubdtype(skeleton_values.dtype, np.integer):
        raise ValueError("skeleton must be integer labels on the image grid")
    if available_values.shape != (*spatial, len(STANDARD13_NAMES)):
        raise ValueError("available must match the image grid and 13 classes")
    if ordering_values.shape != spatial:
        raise ValueError("ordering must match the image grid")
    if not np.isfinite(data).all() or not np.isfinite(centers).all() or not np.isfinite(shapes).all():
        raise ValueError("images, centroids, and transformed shapes must be finite")
    if np.any(labels < 0) or np.any(labels >= len(STANDARD13_NAMES)):
        raise ValueError("segmentation labels must index the standard 13-class profile")
    if not isinstance(source_zero_variance_nan, (bool, np.bool_)):
        raise TypeError("source_zero_variance_nan must be boolean")

    class_count = len(STANDARD13_NAMES)
    relations = lesion_compute_relations(labels, class_count, relation_limit)
    direct_relations = relations.relation_membership[relations.relation_map]

    class_memberships = np.zeros((*spatial, class_count), dtype=np.float32)
    lesion_membership = np.zeros(spatial, dtype=np.float32)
    black_hole_membership = np.zeros(spatial, dtype=np.float32)
    outlier_membership = np.full(spatial, np.float32(0.5), dtype=np.float32)
    lesion_relation = np.ones(spatial, dtype=np.float32)
    object_types = STANDARD13_TYPES
    variances = compute_lesion_toads_variances(
        data,
        class_memberships,
        labels,
        centers,
        lesion_membership,
        black_hole_membership,
        lesion_center,
        black_center,
        object_types,
    )
    channel_weights = compute_lesion_toads_channel_weights(
        data,
        class_memberships,
        centers,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        lesion_center,
        black_center,
        object_types,
        outlier_distance=outlier_energy,
        norm_type=norm_type,
    )

    distance_names: Sequence[str] = STANDARD13_NAMES
    ventricle_mask, ventricle_distance = lesion_toads_ventricle_distance(
        labels, distance_names, maximum_ventricle_distance
    )
    gm_mask, gm_distance = lesion_toads_gm_distance(
        labels, distance_names, maximum_gm_distance
    )
    brainstem_mask, brainstem_distance = lesion_toads_brainstem_distance(
        labels, distance_names, maximum_brainstem_distance
    )
    distances = LesionInitialDistanceMaps(
        ventricle_mask,
        ventricle_distance,
        gm_mask,
        gm_distance,
        brainstem_mask,
        brainstem_distance,
    )

    shape_priors = np.ones(class_count, dtype=np.float32)
    smoothing_priors = np.ones(class_count, dtype=np.float32)
    membership = run_lesion_toads_standard13_volume(
        data,
        centers,
        variances.class_standard_deviation,
        lesion_center,
        variances.lesion_standard_deviation,
        channel_weights,
        shapes,
        shape_priors,
        smoothing_priors,
        direct_relations,
        class_memberships,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        lesion_relation,
        ventricle_distance,
        gm_distance,
        brainstem_distance,
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        gm_flag=True,
        spread=spread,
        source_zero_variance_nan=source_zero_variance_nan,
    )
    loop_entry = seed_lesion_toads_initial_loop(
        membership, iterations=iterations
    )
    first_thinning = run_lesion_initial_thinning_pass(
        loop_entry,
        labels,
        object_types,
        shapes,
        skeleton_values,
        available_values,
        ordering_values,
        lesion_class=10,
        volume_ratio=first_thinning_volume_ratio,
        limit=first_thinning_limit,
    )

    return LesionRegisteredInitialMembership(
        relations,
        variances,
        channel_weights,
        distances,
        class_memberships,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        lesion_relation,
        membership,
        loop_entry,
        first_thinning,
        (
            "constructorState:alreadyOnRegisteredRecropGrid",
            "setGMflag:beforeInitialMemberships",
            "computeRelations",
            "computeVariances",
            "computeChannelWeights",
            "computeVentirclsDistance",
            "computeGMDistance:GMflag=true",
            "computeBstemDistance",
            "computeMemberships:GMflag=true",
            "setGMflag:afterInitialMemberships",
            "Eprev=initialEnergy",
            "Eratio=1.0f",
            "initThinningLabels",
            "propagateThinning",
        ),
        (
            "atlas.updateShapeRegistration:transformAndTransformedShapesAreSupplied",
            "InhomogeneityCorrection2:disabledInSelectedMaCRUISEplusProfile",
        ),
    )
