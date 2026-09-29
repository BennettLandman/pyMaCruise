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
