"""Source-order composition of the selected LesionTOADS membership maps.

This API composes the already-referenced smoothing, prior, relation,
distance-modulation, OUTMASK, and normalization kernels for a single
in-place-order volume pass. Inputs include fixed data-term arrays; estimating
those terms is outside this function. The selected Java method also
accumulates an objective value while updating the maps. The underlying
map-sweep result does not expose that objective, so this adapter makes no
claim about energy/objective equivalence.
"""

from __future__ import annotations

import numpy as np

from .lesion_membership import (
    LesionToadsBlackHoleOutmaskSweepResult,
    run_lesion_toads_blackhole_outmask_smoothing_sweep,
)


def run_lesion_toads_membership_map_pass(
    initial_class_memberships: np.ndarray,
    initial_lesion_membership: np.ndarray,
    initial_black_hole_membership: np.ndarray,
    initial_outlier_membership: np.ndarray,
    class_data_terms: np.ndarray,
    lesion_data_terms: np.ndarray,
    black_hole_data_terms: np.ndarray,
    outmask_data_terms: np.ndarray,
    shape_values: np.ndarray,
    smoothing_priors: np.ndarray,
    class_types: tuple[str, ...] | list[str],
    relation_types: tuple[str, ...] | list[str],
    *,
    smoothing: float,
    shape_priors: np.ndarray | None = None,
    class_similarities: np.ndarray | None = None,
    lesion_similarities: np.ndarray | None = None,
    black_hole_similarities: np.ndarray | None = None,
    prior_coefficient: float = 0.0,
    ventricle_distances: np.ndarray | None = None,
    gray_matter_distances: np.ndarray | None = None,
    brainstem_distances: np.ndarray | None = None,
    max_ventricle_distance: int = 0,
    max_gray_matter_distance: int = 0,
    max_brainstem_distance: int = 0,
    gm_flag: bool = False,
    spread: float = 0.0,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionToadsBlackHoleOutmaskSweepResult:
    """Update and normalize the supported membership maps in source order.

    This is intentionally a small named composition boundary around the
    selected kernels. Supported class types and validation follow
    :func:`run_lesion_toads_blackhole_outmask_smoothing_sweep`: OBJECT,
    optional LESION/OPTIMIZED_LESION, exactly one BLACKHOLE, and exactly one
    OUTMASK. It traverses the mutable maps x, then y, then z, updating each
    voxel before the next one so later voxels observe earlier writes.

    Data terms are supplied as fixed arrays. Statistics estimation, source
    method setup, distance-map generation, unsupported object types, and the
    selected Java objective/energy accumulation are not performed here.
    Arrays are interpreted as RAS-indexed volumes; this low-level function
    has no affine metadata and therefore leaves orientation checking to the
    image-loading boundary.
    """
    return run_lesion_toads_blackhole_outmask_smoothing_sweep(
        initial_class_memberships,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        class_data_terms,
        lesion_data_terms,
        black_hole_data_terms,
        outmask_data_terms,
        shape_values,
        smoothing_priors,
        class_types,
        relation_types,
        smoothing=smoothing,
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        black_hole_similarities=black_hole_similarities,
        prior_coefficient=prior_coefficient,
        ventricle_distances=ventricle_distances,
        gray_matter_distances=gray_matter_distances,
        brainstem_distances=brainstem_distances,
        max_ventricle_distance=max_ventricle_distance,
        max_gray_matter_distance=max_gray_matter_distance,
        max_brainstem_distance=max_brainstem_distance,
        gm_flag=gm_flag,
        spread=spread,
        zero=zero,
        infinity=infinity,
    )
