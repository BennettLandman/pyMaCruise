"""Selected scalar membership conversions from LesionTOADS updates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LesionDistanceModulation:
    """Lesion membership and relation after one source distance chain."""

    membership: np.float32
    relation: np.float32


@dataclass(frozen=True)
class LesionBlackHoleDistanceModulation:
    """Auxiliary lesion maps and relation after BLACKHOLE distance rules."""

    lesion_membership: np.float32
    black_hole_membership: np.float32
    relation: np.float32


@dataclass(frozen=True)
class LesionVoxelMembershipNormalization:
    """Normalized class/auxiliary memberships and energy for one voxel."""

    class_memberships: np.ndarray
    lesion_membership: np.float32
    black_hole_membership: np.float32
    outlier_membership: np.float32
    denominator: np.float32
    energy: np.float32


def convert_lesion_toads_class_membership(
    data_energy: float,
    shape_prior: float,
    *,
    relation_type: str,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> np.float32:
    """Invert one class energy and apply its direct/prior relation rule.

    This reproduces the selected `computeMemberships` class branch before
    denominator accumulation and voxelwise membership normalization. Values
    use Java float32 semantics. Energies at or below the strict ``ZERO``
    threshold map to the selected finite ``INF`` sentinel.
    """
    energy = np.float32(data_energy)
    shape = np.float32(shape_prior)
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(energy):
        raise ValueError("data_energy must be finite")
    if not np.isfinite(shape):
        raise ValueError("shape_prior must be finite")
    if not np.isfinite(zero_value) or zero_value < np.float32(0.0):
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < np.float32(0.0):
        raise ValueError("infinity sentinel must be finite and nonnegative")
    if relation_type not in {"direct", "prior"}:
        raise ValueError("relation_type must be 'direct' or 'prior'")

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        membership = (
            np.float32(np.float32(1.0) / energy)
            if energy > zero_value
            else infinity_value
        )
    if relation_type == "prior":
        doubled_shape = np.float32(np.float32(2.0) * shape)
        # VandyCRUISE Numerics.min(a,b) returns b on ties and when b is NaN.
        shape_factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
        membership = np.float32(shape_factor * membership)
    return membership


def modulate_lesion_toads_membership_by_distances(
    lesion_membership: float,
    lesion_relation: float,
    *,
    ventricle_distance: int,
    gray_matter_distance: int,
    brainstem_distance: int,
    max_ventricle_distance: int,
    max_gray_matter_distance: int,
    max_brainstem_distance: int,
    gm_flag: bool,
    spread: float,
) -> LesionDistanceModulation:
    """Apply the selected lesion/lesion-relation distance modulation.

    A valid ventricle distance takes precedence over valid gray-matter
    distance; a valid brainstem distance is then applied independently. The
    selected source substitutes ``spread`` for zero only in the ventricle and
    gray-matter branches. This helper applies scalar distances already
    computed elsewhere; it does not create masks or distance maps.
    """
    membership = np.float32(lesion_membership)
    relation = np.float32(lesion_relation)
    spread_value = np.float32(spread)
    if not np.isfinite(membership) or not np.isfinite(relation):
        raise ValueError("lesion membership and relation must be finite")
    if not np.isfinite(spread_value):
        raise ValueError("spread must be finite")
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")

    distances = (
        ("ventricle_distance", ventricle_distance),
        ("gray_matter_distance", gray_matter_distance),
        ("brainstem_distance", brainstem_distance),
    )
    maxima = (
        ("max_ventricle_distance", max_ventricle_distance),
        ("max_gray_matter_distance", max_gray_matter_distance),
        ("max_brainstem_distance", max_brainstem_distance),
    )
    for name, value in (*distances, *maxima):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer byte value")
        if not -128 <= int(value) <= 127:
            raise ValueError(f"{name} must fit a signed byte")
    for name, value in distances:
        if int(value) < -1:
            raise ValueError(f"{name} must be -1 or a nonnegative voxel distance")

    vent = int(ventricle_distance)
    gm = int(gray_matter_distance)
    brainstem = int(brainstem_distance)
    vent_max = np.float32(max_ventricle_distance)
    gm_max = np.float32(max_gray_matter_distance)
    brainstem_max = np.float32(max_brainstem_distance)

    if vent_max != np.float32(0.0) and vent != -1:
        numerator = spread_value if vent == 0 else np.float32(vent)
        membership = np.float32(np.float32(numerator * membership) / vent_max)
        factor = np.float32(numerator / vent_max)
        relation = np.float32(relation * factor)
    elif gm_flag and gm_max != np.float32(0.0) and gm != -1:
        numerator = spread_value if gm == 0 else np.float32(gm)
        membership = np.float32(np.float32(numerator * membership) / gm_max)
        factor = np.float32(numerator / gm_max)
        relation = np.float32(relation * factor)

    if brainstem_max != np.float32(0.0) and brainstem != -1:
        numerator = np.float32(brainstem)
        membership = np.float32(np.float32(numerator * membership) / brainstem_max)
        factor = np.float32(numerator / brainstem_max)
        relation = np.float32(relation * factor)

    return LesionDistanceModulation(membership=membership, relation=relation)


def modulate_lesion_toads_blackhole_memberships_by_distances(
    lesion_membership: float,
    black_hole_membership: float,
    lesion_relation: float,
    *,
    ventricle_distance: int,
    gray_matter_distance: int,
    brainstem_distance: int,
    max_ventricle_distance: int,
    max_gray_matter_distance: int,
    max_brainstem_distance: int,
    gm_flag: bool,
    spread: float,
) -> LesionBlackHoleDistanceModulation:
    """Apply the selected BLACKHOLE auxiliary-map distance rules.

    VandyCRUISE's ventricle branch scales ``lesionMembership`` and
    ``lesionRelation`` but leaves ``blackHoleMem`` unchanged. Its GM branch
    scales both auxiliary memberships, and a valid brainstem distance scales
    both independently. The lesion/relation calculations reuse the source-
    checked lesion-only kernel so the selected branch precedence and float32
    order stay aligned with P39z.
    """
    black_hole = np.float32(black_hole_membership)
    if not np.isfinite(black_hole):
        raise ValueError("black_hole_membership must be finite")

    lesion_result = modulate_lesion_toads_membership_by_distances(
        lesion_membership,
        lesion_relation,
        ventricle_distance=ventricle_distance,
        gray_matter_distance=gray_matter_distance,
        brainstem_distance=brainstem_distance,
        max_ventricle_distance=max_ventricle_distance,
        max_gray_matter_distance=max_gray_matter_distance,
        max_brainstem_distance=max_brainstem_distance,
        gm_flag=gm_flag,
        spread=spread,
    )

    # Source BLACKHOLE branch is `if (valid ventricle) ... else if (valid GM)`.
    vent_branch = (
        max_ventricle_distance != 0 and ventricle_distance != -1
    )
    gm_branch = (
        not vent_branch
        and gm_flag
        and max_gray_matter_distance != 0
        and gray_matter_distance != -1
    )
    if gm_branch:
        numerator = (
            np.float32(spread)
            if gray_matter_distance == 0
            else np.float32(gray_matter_distance)
        )
        black_hole = np.float32(
            np.float32(numerator * black_hole)
            / np.float32(max_gray_matter_distance)
        )

    if max_brainstem_distance != 0 and brainstem_distance != -1:
        black_hole = np.float32(
            np.float32(np.float32(brainstem_distance) * black_hole)
            / np.float32(max_brainstem_distance)
        )

    return LesionBlackHoleDistanceModulation(
        lesion_membership=lesion_result.membership,
        black_hole_membership=black_hole,
        relation=lesion_result.relation,
    )


def normalize_lesion_toads_voxel_memberships(
    class_memberships: np.ndarray,
    class_energies: np.ndarray,
    class_types: tuple[str, ...] | list[str],
    *,
    lesion_membership: float,
    black_hole_membership: float,
    outlier_membership: float,
    lesion_energy: float,
    black_hole_energy: float,
) -> LesionVoxelMembershipNormalization:
    """Run the selected final denominator/normalization loop for one voxel.

    Inputs are already-assigned class and auxiliary memberships plus their
    local energies. This deliberately preserves source traversal behavior,
    including adding lesion auxiliaries to the denominator for each matching
    class and dividing shared auxiliary maps once per matching class during
    normalization.
    """
    memberships = np.asarray(class_memberships, dtype=np.float32).copy()
    energies = np.asarray(class_energies, dtype=np.float32)
    types = tuple(class_types)
    if memberships.ndim != 1 or energies.ndim != 1:
        raise ValueError("class_memberships and class_energies must be vectors")
    if len(memberships) == 0 or len(energies) != len(memberships):
        raise ValueError("class membership and energy vectors must have equal nonzero length")
    if len(types) != len(memberships):
        raise ValueError("class_types must have one entry per class")
    supported = {"OBJECT", "LESION", "OPTIMIZED_LESION", "BLACKHOLE", "OUTLIER", "OUTMASK", "MASK"}
    if any(not isinstance(kind, str) or kind not in supported for kind in types):
        raise ValueError("class_types contains an unsupported LesionTOADS type")

    lesion = np.float32(lesion_membership)
    black_hole = np.float32(black_hole_membership)
    outlier = np.float32(outlier_membership)
    lesion_term = np.float32(lesion_energy)
    black_hole_term = np.float32(black_hole_energy)

    denominator = np.float32(0.0)
    for membership, kind in zip(memberships, types):
        if kind in {"LESION", "OPTIMIZED_LESION"}:
            denominator = np.float32(denominator + lesion)
        elif kind == "BLACKHOLE":
            auxiliary_sum = np.float32(lesion + black_hole)
            denominator = np.float32(denominator + auxiliary_sum)
        denominator = np.float32(denominator + membership)

    local_energy = np.float32(0.0)
    if denominator > np.float32(0.0):
        for index, kind in enumerate(types):
            memberships[index] = np.float32(memberships[index] / denominator)
            if kind in {"LESION", "OPTIMIZED_LESION"}:
                lesion = np.float32(lesion / denominator)
            elif kind == "BLACKHOLE":
                lesion = np.float32(lesion / denominator)
                black_hole = np.float32(black_hole / denominator)
            elif kind in {"OUTLIER", "OUTMASK"}:
                outlier = np.float32(outlier / denominator)

            energy_membership = np.float32(energies[index] * memberships[index])
            energy_membership = np.float32(energy_membership * memberships[index])
            local_energy = np.float32(local_energy + energy_membership)

    lesion_energy_value = np.float32(lesion_term * lesion)
    lesion_energy_value = np.float32(lesion_energy_value * lesion)
    black_hole_energy_value = np.float32(black_hole_term * black_hole)
    black_hole_energy_value = np.float32(black_hole_energy_value * black_hole)
    auxiliary_energy = np.float32(lesion_energy_value + black_hole_energy_value)
    local_energy = np.float32(local_energy + auxiliary_energy)

    return LesionVoxelMembershipNormalization(
        class_memberships=memberships,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        outlier_membership=outlier,
        denominator=denominator,
        energy=local_energy,
    )
