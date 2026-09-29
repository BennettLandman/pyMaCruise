"""Selected scalar membership conversions from LesionTOADS updates."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Iterator

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


@dataclass(frozen=True)
class LesionToadsSmoothingEnergyUpdate:
    """One class/auxiliary energy update from the selected six-face stencil."""

    class_membership_data_term: np.float32
    lesion_membership_data_term: np.float32
    black_hole_membership_data_term: np.float32
    outmask_membership_data_term: np.float32
    class_energy: np.float32
    lesion_energy: np.float32
    black_hole_energy: np.float32
    other_class_neighbor_squares: np.float32
    lesion_neighbor_squares: np.float32
    black_hole_neighbor_squares: np.float32
    current_class_neighbor_squares: np.float32


@dataclass(frozen=True)
class LesionToadsPriorEnergyUpdate:
    """One voxel's pairwise, lesion, and black-hole prior contributions."""

    class_prior: np.float32
    lesion_prior: np.float32
    black_hole_prior: np.float32
    class_membership_data_term: np.float32
    lesion_membership_data_term: np.float32
    black_hole_membership_data_term: np.float32
    outmask_membership_data_term: np.float32
    class_energy: np.float32
    lesion_energy: np.float32
    black_hole_energy: np.float32


@dataclass(frozen=True)
class LesionToadsOutmaskMembership:
    """OUTMASK class and outlier-map values after direct/prior mixing."""

    class_membership: np.float32
    outlier_membership: np.float32


@dataclass(frozen=True)
class LesionToadsSixFaceNeighborhood:
    """One interior voxel and its source-ordered six-face membership values."""

    coordinate: tuple[int, int, int]
    class_memberships: np.ndarray
    lesion_memberships: np.ndarray
    black_hole_memberships: np.ndarray


@dataclass(frozen=True)
class LesionToadsVolumeSmoothingUpdate:
    """Per-class smoothing-energy updates at one streamed volume location."""

    coordinate: tuple[int, int, int]
    class_updates: tuple[LesionToadsSmoothingEnergyUpdate, ...]


@dataclass(frozen=True)
class LesionToadsLesionSweepResult:
    """Class and auxiliary lesion memberships after the bounded volume sweep."""

    class_memberships: np.ndarray
    lesion_membership: np.ndarray


@dataclass(frozen=True)
class LesionToadsBlackHoleSweepResult:
    """Class, lesion, and black-hole maps after a bounded volume sweep."""

    class_memberships: np.ndarray
    lesion_membership: np.ndarray
    black_hole_membership: np.ndarray


@dataclass(frozen=True)
class LesionToadsOutmaskSweepResult:
    """Class and outlier maps after the bounded OUTMASK volume sweep."""

    class_memberships: np.ndarray
    outlier_membership: np.ndarray


@dataclass(frozen=True)
class LesionToadsLesionOutmaskSweepResult:
    """Class, lesion, and outlier maps after a bounded combined sweep."""

    class_memberships: np.ndarray
    lesion_membership: np.ndarray
    outlier_membership: np.ndarray


@dataclass(frozen=True)
class LesionToadsBlackHoleOutmaskSweepResult:
    """Class, lesion, black-hole, and outlier maps after a combined sweep."""

    class_memberships: np.ndarray
    lesion_membership: np.ndarray
    black_hole_membership: np.ndarray
    outlier_membership: np.ndarray


def iter_lesion_toads_six_face_neighborhoods(
    class_memberships: np.ndarray,
    lesion_memberships: np.ndarray,
    black_hole_memberships: np.ndarray,
) -> Iterator[LesionToadsSixFaceNeighborhood]:
    """Yield source-ordered six-face neighborhoods for interior voxels.

    The selected VandyCRUISE update loops x, then y, then z from 1 through
    size-2 and reads faces in X+, X-, Y+, Y-, Z+, Z- order. Boundary voxels
    are not yielded, matching the source loop's unchanged-border behavior.
    Inputs are treated as read-only; neighborhoods are streamed to avoid a
    large temporary stencil tensor on full atlas-sized volumes.
    """
    classes = np.asarray(class_memberships, dtype=np.float32)
    lesion = np.asarray(lesion_memberships, dtype=np.float32)
    black_hole = np.asarray(black_hole_memberships, dtype=np.float32)
    if classes.ndim != 4 or any(size < 3 for size in classes.shape[:3]) or classes.shape[3] == 0:
        raise ValueError("class_memberships must have shape (x, y, z, classes), with each spatial size >= 3")
    if lesion.shape != classes.shape[:3] or black_hole.shape != classes.shape[:3]:
        raise ValueError("lesion and black-hole memberships must match the class spatial shape")

    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for x in range(1, classes.shape[0] - 1):
        for y in range(1, classes.shape[1] - 1):
            for z in range(1, classes.shape[2] - 1):
                faces = tuple((x + dx, y + dy, z + dz) for dx, dy, dz in offsets)
                yield LesionToadsSixFaceNeighborhood(
                    coordinate=(x, y, z),
                    class_memberships=np.stack([classes[point] for point in faces]),
                    lesion_memberships=np.asarray([lesion[point] for point in faces], dtype=np.float32),
                    black_hole_memberships=np.asarray([black_hole[point] for point in faces], dtype=np.float32),
                )


def iter_lesion_toads_volume_smoothing_updates(
    class_memberships: np.ndarray,
    lesion_memberships: np.ndarray,
    black_hole_memberships: np.ndarray,
    class_membership_data_terms: np.ndarray,
    lesion_membership_data_terms: np.ndarray,
    black_hole_membership_data_terms: np.ndarray,
    outmask_membership_data_terms: np.ndarray,
    class_energies: np.ndarray,
    lesion_energies: np.ndarray,
    black_hole_energies: np.ndarray,
    class_types: tuple[str, ...] | list[str],
    smoothing_priors: np.ndarray,
    *,
    smoothing: float,
) -> Iterator[LesionToadsVolumeSmoothingUpdate]:
    """Stream the selected smoothing kernel over fixed synthetic/source maps.

    This composes neighborhood extraction with the already validated scalar
    update. Term/energy maps are one value per voxel and class. It intentionally
    does not write memberships, apply atlas priors, invert terms, or normalize;
    a later full update must reproduce Java's in-place scan behavior.
    """
    memberships = np.asarray(class_memberships, dtype=np.float32)
    if memberships.ndim != 4 or memberships.shape[3] == 0:
        raise ValueError("class_memberships must have shape (x, y, z, classes)")
    spatial_shape = memberships.shape[:3]
    classes = memberships.shape[3]
    types = tuple(class_types)
    if len(types) != classes:
        raise ValueError("class_types must have one entry per class")
    term_names = (
        "class_membership_data_terms",
        "lesion_membership_data_terms",
        "black_hole_membership_data_terms",
        "outmask_membership_data_terms",
        "class_energies",
        "lesion_energies",
        "black_hole_energies",
    )
    terms = tuple(
        np.asarray(value, dtype=np.float32)
        for value in (
            class_membership_data_terms,
            lesion_membership_data_terms,
            black_hole_membership_data_terms,
            outmask_membership_data_terms,
            class_energies,
            lesion_energies,
            black_hole_energies,
        )
    )
    for name, values in zip(term_names, terms, strict=True):
        if values.shape != (*spatial_shape, classes):
            raise ValueError(f"{name} must match class_memberships shape")
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (classes,):
        raise ValueError("smoothing_priors must have one value per class")
    if not np.isfinite(priors).all() or not np.isfinite(np.float32(smoothing)):
        raise ValueError("smoothing parameters must be finite")
    lesion = np.asarray(lesion_memberships, dtype=np.float32)
    black_hole = np.asarray(black_hole_memberships, dtype=np.float32)
    if lesion.shape != spatial_shape or black_hole.shape != spatial_shape:
        raise ValueError("lesion and black-hole memberships must match the class spatial shape")

    for neighborhood in iter_lesion_toads_six_face_neighborhoods(memberships, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        updates = tuple(
            update_lesion_toads_smoothing_energies(
                neighborhood.class_memberships,
                neighborhood.lesion_memberships,
                neighborhood.black_hole_memberships,
                class_index=k,
                class_type=types[k],
                smoothing_prior=float(priors[k]),
                smoothing=smoothing,
                class_energy=float(terms[4][x, y, z, k]),
                lesion_energy=float(terms[5][x, y, z, k]),
                black_hole_energy=float(terms[6][x, y, z, k]),
                class_membership_data_term=float(terms[0][x, y, z, k]),
                lesion_membership_data_term=float(terms[1][x, y, z, k]),
                black_hole_membership_data_term=float(terms[2][x, y, z, k]),
                outmask_membership_data_term=float(terms[3][x, y, z, k]),
            )
            for k in range(classes)
        )
        yield LesionToadsVolumeSmoothingUpdate(neighborhood.coordinate, updates)


def run_lesion_toads_object_membership_sweep(
    initial_class_memberships: np.ndarray,
    class_data_terms: np.ndarray,
    smoothing_priors: np.ndarray,
    *,
    smoothing: float,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> np.ndarray:
    """Run the direct, class-only OBJECT subset of the source volume sweep.

    All classes are assumed to have positive shape support and direct
    relations; per-voxel terms contain only the supplied data term plus the
    six-face smoothing addition. Atlas priors, masks, auxiliary classes, and
    lesion-specific distance handling are deliberately excluded. The method
    mutates a private copy in Java's x/y/z and class order, normalizes each
    voxel before advancing, and therefore lets later voxels observe prior
    updates exactly as the selected in-place scan does.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    if data.shape != current.shape:
        raise ValueError("class_data_terms must match initial_class_memberships shape")
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (current.shape[3],) or not np.isfinite(priors).all() or np.any(priors < 0):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not np.isfinite(current).all() or not np.isfinite(data).all():
        raise ValueError("membership and data-term volumes must be finite")
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= np.float32(0.0):
        raise ValueError("smoothing must be finite and positive")

    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    lesion = np.zeros(current.shape[:3], dtype=np.float32)
    black_hole = np.zeros(current.shape[:3], dtype=np.float32)
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(current.shape[3]):
            update = update_lesion_toads_smoothing_energies(
                neighborhood.class_memberships,
                neighborhood.lesion_memberships,
                neighborhood.black_hole_memberships,
                class_index=k,
                class_type="OBJECT",
                smoothing_prior=float(priors[k]),
                smoothing=float(smoothing_value),
                class_energy=0.0,
                lesion_energy=0.0,
                black_hole_energy=0.0,
                class_membership_data_term=float(data[x, y, z, k]),
                lesion_membership_data_term=0.0,
                black_hole_membership_data_term=0.0,
                outmask_membership_data_term=0.0,
            )
            current[x, y, z, k] = convert_lesion_toads_class_membership(
                float(update.class_membership_data_term),
                1.0,
                relation_type="direct",
                zero=float(zero_value),
                infinity=float(infinity_value),
            )
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(current.shape[3]):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
    return current


def run_lesion_toads_object_prior_membership_sweep(
    initial_class_memberships: np.ndarray,
    class_data_terms: np.ndarray,
    shape_values: np.ndarray,
    shape_priors: np.ndarray,
    class_similarities: np.ndarray,
    smoothing_priors: np.ndarray,
    relation_types: tuple[str, ...] | list[str],
    *,
    smoothing: float,
    prior_coefficient: float,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> np.ndarray:
    """Sweep ordinary OBJECT classes with source shape-prior energy terms.

    This keeps the P39zo restrictions (direct object classes and no lesion or
    auxiliary maps) while adding the selected pairwise atlas prior. Shape
    values and class relations are caller-supplied per voxel/class; the source
    class-prior coefficients are a fixed class-by-class matrix. Memberships
    are updated and normalized in place on a private copy in x/y/z order.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    relations = tuple(relation_types)
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    priors = np.asarray(shape_priors, dtype=np.float32)
    similarities = np.asarray(class_similarities, dtype=np.float32)
    smoothing_values = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (classes,) or smoothing_values.shape != (classes,):
        raise ValueError("shape_priors and smoothing_priors must have one value per class")
    if similarities.shape != (classes, classes):
        raise ValueError("class_similarities must have shape (classes, classes)")
    if not all(np.isfinite(array).all() for array in (current, data, shapes, priors, similarities, smoothing_values)):
        raise ValueError("membership, term, and prior arrays must be finite")
    if np.any(priors < 0) or np.any(smoothing_values < 0):
        raise ValueError("shape and smoothing priors must be nonnegative")
    smoothing_value = np.float32(smoothing)
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    if not np.isfinite(coefficient) or coefficient < 0:
        raise ValueError("prior_coefficient must be finite and nonnegative")

    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")
    lesion = np.zeros(spatial_shape, dtype=np.float32)
    black_hole = np.zeros(spatial_shape, dtype=np.float32)
    zero_similarity = np.zeros(classes, dtype=np.float32)

    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            shape = shapes[x, y, z]
            if shape[k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
                continue
            update = update_lesion_toads_smoothing_energies(
                neighborhood.class_memberships,
                neighborhood.lesion_memberships,
                neighborhood.black_hole_memberships,
                class_index=k,
                class_type="OBJECT",
                smoothing_prior=float(smoothing_values[k]),
                smoothing=float(smoothing_value),
                class_energy=0.0,
                lesion_energy=0.0,
                black_hole_energy=0.0,
                class_membership_data_term=float(data[x, y, z, k]),
                lesion_membership_data_term=0.0,
                black_hole_membership_data_term=0.0,
                outmask_membership_data_term=0.0,
            )
            class_term = np.float32(update.class_membership_data_term)
            if coefficient > np.float32(0.0):
                prior_update = update_lesion_toads_prior_energies(
                    priors,
                    similarities[k],
                    zero_similarity,
                    zero_similarity,
                    shape,
                    class_index=k,
                    class_type="OBJECT",
                    prior_coefficient=float(coefficient),
                    class_membership_data_term=float(class_term),
                    lesion_membership_data_term=0.0,
                    black_hole_membership_data_term=0.0,
                    outmask_membership_data_term=0.0,
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    has_lesion=False,
                    lesion_class_index=0,
                    lesion_class_type="LESION",
                )
                class_term = np.float32(prior_update.class_membership_data_term)
            current[x, y, z, k] = convert_lesion_toads_class_membership(
                float(class_term),
                float(shape[k]),
                relation_type=relations[k],
                zero=float(zero_value),
                infinity=float(infinity_value),
            )
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
    return current


def run_lesion_toads_lesion_smoothing_sweep(
    initial_class_memberships: np.ndarray,
    initial_lesion_membership: np.ndarray,
    class_data_terms: np.ndarray,
    lesion_data_terms: np.ndarray,
    shape_values: np.ndarray,
    smoothing_priors: np.ndarray,
    class_types: tuple[str, ...] | list[str],
    relation_types: tuple[str, ...] | list[str],
    *,
    smoothing: float,
    shape_priors: np.ndarray | None = None,
    class_similarities: np.ndarray | None = None,
    lesion_similarities: np.ndarray | None = None,
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
) -> LesionToadsLesionSweepResult:
    """Sweep OBJECT classes plus one LESION class and its auxiliary map.

    The selected LESION branch receives a separate lesion data term and
    updates `lesionMembership` from the same smoothing stencil. The method
    preserves source shape gating, class relation factors, in-place x/y/z
    neighbor reads, denominator accumulation (auxiliary lesion before the
    lesion class), and per-voxel normalization. Optional class/lesion atlas
    priors reuse the selected voxel kernel; distance maps and other auxiliary
    or class types remain outside this helper. Optional distance maps apply
    after lesion conversion and before the shared denominator.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("LESION") != 1 or any(
        kind not in {"OBJECT", "LESION"} for kind in types
    ):
        raise ValueError("class_types must contain ordinary OBJECT classes and exactly one LESION class")
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    if lesion.shape != spatial_shape or lesion_data.shape != spatial_shape:
        raise ValueError("lesion membership and data-term maps must match the class spatial shape")
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (classes,) or not np.isfinite(priors).all() or np.any(priors < 0):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not all(np.isfinite(array).all() for array in (current, lesion, data, lesion_data, shapes)):
        raise ValueError("membership, data-term, and shape arrays must be finite")
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient) or coefficient < np.float32(0.0):
        raise ValueError("prior_coefficient must be finite and nonnegative")
    atlas_priors = class_similarity = lesion_similarity = None
    if coefficient > np.float32(0.0):
        if shape_priors is None or class_similarities is None or lesion_similarities is None:
            raise ValueError("shape priors and class/lesion similarities are required for a positive prior_coefficient")
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        class_similarity = np.asarray(class_similarities, dtype=np.float32)
        lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or lesion_similarity.shape != (classes,):
            raise ValueError("shape_priors and lesion_similarities must have one value per class")
        if class_similarity.shape != (classes, classes):
            raise ValueError("class_similarities must have shape (classes, classes)")
        if not all(np.isfinite(array).all() for array in (atlas_priors, class_similarity, lesion_similarity)):
            raise ValueError("atlas prior and similarity inputs must be finite")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    maxima = (
        ("max_ventricle_distance", max_ventricle_distance),
        ("max_gray_matter_distance", max_gray_matter_distance),
        ("max_brainstem_distance", max_brainstem_distance),
    )
    for name, value in maxima:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not -128 <= int(value) <= 127:
            raise ValueError(f"{name} must fit a signed byte")
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")
    spread_value = np.float32(spread)
    if not np.isfinite(spread_value):
        raise ValueError("spread must be finite")
    distance_inputs = (
        ("ventricle_distances", ventricle_distances, max_ventricle_distance),
        ("gray_matter_distances", gray_matter_distances, max_gray_matter_distance),
        ("brainstem_distances", brainstem_distances, max_brainstem_distance),
    )
    distance_maps: dict[str, np.ndarray] = {}
    for name, source, maximum in distance_inputs:
        if source is None:
            if maximum != 0:
                raise ValueError(f"{name} is required when its maximum distance is nonzero")
            distance_maps[name] = np.full(spatial_shape, -1, dtype=np.int8)
            continue
        raw = np.asarray(source)
        if not np.issubdtype(raw.dtype, np.integer) or raw.shape != spatial_shape:
            raise ValueError(f"{name} must be an integer map matching the class spatial shape")
        if np.any(raw < -1) or np.any(raw > 127):
            raise ValueError(f"{name} values must fit the source -1..127 distance range")
        distance_maps[name] = np.asarray(raw, dtype=np.int8)
    use_distances = any(value != 0 for _, value in maxima)

    black_hole = np.zeros(spatial_shape, dtype=np.float32)
    zero_similarity = np.zeros(classes, dtype=np.float32)
    lesion_class = types.index("LESION")
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            if shapes[x, y, z, k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
                if k == lesion_class:
                    lesion[x, y, z] = np.float32(0.0)
            else:
                update = update_lesion_toads_smoothing_energies(
                    neighborhood.class_memberships,
                    neighborhood.lesion_memberships,
                    neighborhood.black_hole_memberships,
                    class_index=k,
                    class_type=types[k],
                    smoothing_prior=float(priors[k]),
                    smoothing=float(smoothing_value),
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    class_membership_data_term=float(data[x, y, z, k]),
                    lesion_membership_data_term=(float(lesion_data[x, y, z]) if k == lesion_class else 0.0),
                    black_hole_membership_data_term=0.0,
                    outmask_membership_data_term=0.0,
                )
                class_term = np.float32(update.class_membership_data_term)
                lesion_term = np.float32(update.lesion_membership_data_term)
                if coefficient > np.float32(0.0):
                    prior_update = update_lesion_toads_prior_energies(
                        atlas_priors,
                        class_similarity[k],
                        lesion_similarity,
                        zero_similarity,
                        shapes[x, y, z],
                        class_index=k,
                        class_type=types[k],
                        prior_coefficient=float(coefficient),
                        class_membership_data_term=float(class_term),
                        lesion_membership_data_term=float(lesion_term),
                        black_hole_membership_data_term=0.0,
                        outmask_membership_data_term=0.0,
                        class_energy=0.0,
                        lesion_energy=0.0,
                        black_hole_energy=0.0,
                        has_lesion=True,
                        lesion_class_index=lesion_class,
                        lesion_class_type="LESION",
                    )
                    class_term = np.float32(prior_update.class_membership_data_term)
                    lesion_term = np.float32(prior_update.lesion_membership_data_term)
                current[x, y, z, k] = convert_lesion_toads_class_membership(
                    float(class_term),
                    float(shapes[x, y, z, k]),
                    relation_type=relations[k],
                    zero=float(zero_value),
                    infinity=float(infinity_value),
                )
                if k == lesion_class:
                    lesion[x, y, z] = convert_lesion_toads_class_membership(
                        float(lesion_term),
                        float(shapes[x, y, z, k]),
                        relation_type=relations[k],
                        zero=float(zero_value),
                        infinity=float(infinity_value),
                    )
                    if use_distances:
                        shape = np.float32(shapes[x, y, z, k])
                        lesion_relation = (
                            np.float32(1.0)
                            if relations[k] == "direct"
                            else np.float32(1.0)
                            if np.float32(1.0) < np.float32(np.float32(2.0) * shape)
                            else np.float32(np.float32(2.0) * shape)
                        )
                        modulated = modulate_lesion_toads_membership_by_distances(
                            float(lesion[x, y, z]),
                            float(lesion_relation),
                            ventricle_distance=int(distance_maps["ventricle_distances"][x, y, z]),
                            gray_matter_distance=int(distance_maps["gray_matter_distances"][x, y, z]),
                            brainstem_distance=int(distance_maps["brainstem_distances"][x, y, z]),
                            max_ventricle_distance=int(max_ventricle_distance),
                            max_gray_matter_distance=int(max_gray_matter_distance),
                            max_brainstem_distance=int(max_brainstem_distance),
                            gm_flag=gm_flag,
                            spread=float(spread_value),
                        )
                        lesion[x, y, z] = modulated.membership
            if k == lesion_class:
                denominator = np.float32(denominator + lesion[x, y, z])
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if k == lesion_class:
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
    return LesionToadsLesionSweepResult(current, lesion)


def run_lesion_toads_blackhole_smoothing_sweep(
    initial_class_memberships: np.ndarray,
    initial_lesion_membership: np.ndarray,
    initial_black_hole_membership: np.ndarray,
    class_data_terms: np.ndarray,
    lesion_data_terms: np.ndarray,
    black_hole_data_terms: np.ndarray,
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
) -> LesionToadsBlackHoleSweepResult:
    """Sweep ordinary OBJECT classes plus one BLACKHOLE auxiliary branch.

    This composes the selected BLACKHOLE six-face smoothing equations with
    separate lesion and black-hole data terms, class/auxiliary inversion,
    source denominator order, and per-voxel normalization. Optional
    class/lesion/black-hole atlas priors reuse the selected voxel kernel;
    optional black-hole-specific distance modulation uses the source's
    ventricle/GM/brainstem asymmetry.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    black_hole = np.asarray(initial_black_hole_membership, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    black_hole_data = np.asarray(black_hole_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("BLACKHOLE") != 1 or any(
        kind not in {"OBJECT", "BLACKHOLE"} for kind in types
    ):
        raise ValueError("class_types must contain ordinary OBJECT classes and exactly one BLACKHOLE class")
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    if any(array.shape != spatial_shape for array in (lesion, black_hole, lesion_data, black_hole_data)):
        raise ValueError("lesion/black-hole membership and data maps must match the class spatial shape")
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (classes,) or not np.isfinite(priors).all() or np.any(priors < 0):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not all(np.isfinite(array).all() for array in (current, lesion, black_hole, data, lesion_data, black_hole_data, shapes)):
        raise ValueError("membership, data-term, and shape arrays must be finite")
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient) or coefficient < np.float32(0.0):
        raise ValueError("prior_coefficient must be finite and nonnegative")
    atlas_priors = class_similarity = lesion_similarity = black_hole_similarity = None
    if coefficient > np.float32(0.0):
        if any(value is None for value in (shape_priors, class_similarities, lesion_similarities, black_hole_similarities)):
            raise ValueError("shape priors and class/lesion/black-hole similarities are required for a positive prior_coefficient")
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        class_similarity = np.asarray(class_similarities, dtype=np.float32)
        lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
        black_hole_similarity = np.asarray(black_hole_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or lesion_similarity.shape != (classes,) or black_hole_similarity.shape != (classes,):
            raise ValueError("shape priors and lesion/black-hole similarities must have one value per class")
        if class_similarity.shape != (classes, classes):
            raise ValueError("class_similarities must have shape (classes, classes)")
        if not all(np.isfinite(array).all() for array in (atlas_priors, class_similarity, lesion_similarity, black_hole_similarity)):
            raise ValueError("atlas prior and similarity inputs must be finite")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    maxima = (
        ("max_ventricle_distance", max_ventricle_distance),
        ("max_gray_matter_distance", max_gray_matter_distance),
        ("max_brainstem_distance", max_brainstem_distance),
    )
    for name, value in maxima:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not -128 <= int(value) <= 127:
            raise ValueError(f"{name} must fit a signed byte")
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")
    spread_value = np.float32(spread)
    if not np.isfinite(spread_value):
        raise ValueError("spread must be finite")
    distance_inputs = (
        ("ventricle_distances", ventricle_distances, max_ventricle_distance),
        ("gray_matter_distances", gray_matter_distances, max_gray_matter_distance),
        ("brainstem_distances", brainstem_distances, max_brainstem_distance),
    )
    distance_maps: dict[str, np.ndarray] = {}
    for name, source, maximum in distance_inputs:
        if source is None:
            if maximum != 0:
                raise ValueError(f"{name} is required when its maximum distance is nonzero")
            distance_maps[name] = np.full(spatial_shape, -1, dtype=np.int8)
            continue
        raw = np.asarray(source)
        if not np.issubdtype(raw.dtype, np.integer) or raw.shape != spatial_shape:
            raise ValueError(f"{name} must be an integer map matching the class spatial shape")
        if np.any(raw < -1) or np.any(raw > 127):
            raise ValueError(f"{name} values must fit the source -1..127 distance range")
        distance_maps[name] = np.asarray(raw, dtype=np.int8)
    use_distances = any(value != 0 for _, value in maxima)

    black_hole_class = types.index("BLACKHOLE")
    zero_similarity = np.zeros(classes, dtype=np.float32)
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            if shapes[x, y, z, k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
                if k == black_hole_class:
                    lesion[x, y, z] = np.float32(0.0)
                    black_hole[x, y, z] = np.float32(0.0)
            else:
                update = update_lesion_toads_smoothing_energies(
                    neighborhood.class_memberships,
                    neighborhood.lesion_memberships,
                    neighborhood.black_hole_memberships,
                    class_index=k,
                    class_type=types[k],
                    smoothing_prior=float(priors[k]),
                    smoothing=float(smoothing_value),
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    class_membership_data_term=float(data[x, y, z, k]),
                    lesion_membership_data_term=(float(lesion_data[x, y, z]) if k == black_hole_class else 0.0),
                    black_hole_membership_data_term=(float(black_hole_data[x, y, z]) if k == black_hole_class else 0.0),
                    outmask_membership_data_term=0.0,
                )
                class_term = np.float32(update.class_membership_data_term)
                lesion_term = np.float32(update.lesion_membership_data_term)
                black_hole_term = np.float32(update.black_hole_membership_data_term)
                if coefficient > np.float32(0.0):
                    prior_update = update_lesion_toads_prior_energies(
                        atlas_priors,
                        class_similarity[k],
                        lesion_similarity,
                        black_hole_similarity,
                        shapes[x, y, z],
                        class_index=k,
                        class_type=types[k],
                        prior_coefficient=float(coefficient),
                        class_membership_data_term=float(class_term),
                        lesion_membership_data_term=float(lesion_term),
                        black_hole_membership_data_term=float(black_hole_term),
                        outmask_membership_data_term=0.0,
                        class_energy=0.0,
                        lesion_energy=0.0,
                        black_hole_energy=0.0,
                        has_lesion=True,
                        lesion_class_index=black_hole_class,
                        lesion_class_type="BLACKHOLE",
                    )
                    class_term = np.float32(prior_update.class_membership_data_term)
                    lesion_term = np.float32(prior_update.lesion_membership_data_term)
                    black_hole_term = np.float32(prior_update.black_hole_membership_data_term)
                current[x, y, z, k] = convert_lesion_toads_class_membership(
                    float(class_term),
                    float(shapes[x, y, z, k]),
                    relation_type=relations[k],
                    zero=float(zero_value),
                    infinity=float(infinity_value),
                )
                if k == black_hole_class:
                    lesion[x, y, z] = convert_lesion_toads_class_membership(
                        float(lesion_term),
                        float(shapes[x, y, z, k]),
                        relation_type=relations[k],
                        zero=float(zero_value),
                        infinity=float(infinity_value),
                    )
                    black_hole[x, y, z] = convert_lesion_toads_class_membership(
                        float(black_hole_term),
                        float(shapes[x, y, z, k]),
                        relation_type=relations[k],
                        zero=float(zero_value),
                        infinity=float(infinity_value),
                    )
                    if use_distances:
                        shape = np.float32(shapes[x, y, z, k])
                        lesion_relation = (
                            np.float32(1.0)
                            if relations[k] == "direct"
                            else np.float32(1.0)
                            if np.float32(1.0) < np.float32(np.float32(2.0) * shape)
                            else np.float32(np.float32(2.0) * shape)
                        )
                        modulated = modulate_lesion_toads_blackhole_memberships_by_distances(
                            float(lesion[x, y, z]),
                            float(black_hole[x, y, z]),
                            float(lesion_relation),
                            ventricle_distance=int(distance_maps["ventricle_distances"][x, y, z]),
                            gray_matter_distance=int(distance_maps["gray_matter_distances"][x, y, z]),
                            brainstem_distance=int(distance_maps["brainstem_distances"][x, y, z]),
                            max_ventricle_distance=int(max_ventricle_distance),
                            max_gray_matter_distance=int(max_gray_matter_distance),
                            max_brainstem_distance=int(max_brainstem_distance),
                            gm_flag=gm_flag,
                            spread=float(spread_value),
                        )
                        lesion[x, y, z] = modulated.lesion_membership
                        black_hole[x, y, z] = modulated.black_hole_membership
            if k == black_hole_class:
                auxiliary_sum = np.float32(lesion[x, y, z] + black_hole[x, y, z])
                denominator = np.float32(denominator + auxiliary_sum)
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if k == black_hole_class:
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                    black_hole[x, y, z] = np.float32(black_hole[x, y, z] / denominator)
    return LesionToadsBlackHoleSweepResult(current, lesion, black_hole)


def run_lesion_toads_outmask_smoothing_sweep(
    initial_class_memberships: np.ndarray,
    initial_outlier_membership: np.ndarray,
    class_data_terms: np.ndarray,
    outmask_data_terms: np.ndarray,
    shape_values: np.ndarray,
    smoothing_priors: np.ndarray,
    class_types: tuple[str, ...] | list[str],
    relation_types: tuple[str, ...] | list[str],
    *,
    smoothing: float,
    shape_priors: np.ndarray | None = None,
    class_similarities: np.ndarray | None = None,
    prior_coefficient: float = 0.0,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionToadsOutmaskSweepResult:
    """Sweep ordinary OBJECT classes plus one OUTMASK/outlier map.

    This composes selected six-face smoothing with inversion of the OUTMASK
    and outlier terms, direct/prior blending, in-place x/y/z neighborhood
    reads, and source denominator normalization. Optional atlas similarities
    add selected class-prior terms to both the class and OUTMASK energies.
    Other auxiliary/class branches remain separate.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    outmask_data = np.asarray(outmask_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("OUTMASK") != 1 or any(
        kind not in {"OBJECT", "OUTMASK"} for kind in types
    ):
        raise ValueError("class_types must contain ordinary OBJECT classes and exactly one OUTMASK class")
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    if outlier.shape != spatial_shape or outmask_data.shape != spatial_shape:
        raise ValueError("outlier membership and OUTMASK data-term maps must match the class spatial shape")
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    if priors.shape != (classes,) or not np.isfinite(priors).all() or np.any(priors < 0):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not all(np.isfinite(array).all() for array in (current, outlier, data, outmask_data, shapes)):
        raise ValueError("membership, data-term, and shape arrays must be finite")
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient) or coefficient < np.float32(0.0):
        raise ValueError("prior_coefficient must be finite and nonnegative")
    if shape_priors is None and class_similarities is None:
        if coefficient != np.float32(0.0):
            raise ValueError("shape_priors and class_similarities are required for a nonzero prior coefficient")
        atlas_priors = np.zeros(classes, dtype=np.float32)
        similarity = np.zeros((classes, classes), dtype=np.float32)
    elif shape_priors is None or class_similarities is None:
        raise ValueError("shape_priors and class_similarities must be provided together")
    else:
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        similarity = np.asarray(class_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or similarity.shape != (classes, classes):
            raise ValueError("shape_priors and class_similarities must have shapes (classes,) and (classes, classes)")
        if not np.isfinite(atlas_priors).all() or not np.isfinite(similarity).all():
            raise ValueError("shape priors and class similarities must be finite")
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    zero_map = np.zeros(spatial_shape, dtype=np.float32)
    zero_similarities = np.zeros(classes, dtype=np.float32)
    outmask_class = types.index("OUTMASK")
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, zero_map, zero_map):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            if shapes[x, y, z, k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
            else:
                update = update_lesion_toads_smoothing_energies(
                    neighborhood.class_memberships,
                    neighborhood.lesion_memberships,
                    neighborhood.black_hole_memberships,
                    class_index=k,
                    class_type=types[k],
                    smoothing_prior=float(priors[k]),
                    smoothing=float(smoothing_value),
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    class_membership_data_term=float(data[x, y, z, k]),
                    lesion_membership_data_term=0.0,
                    black_hole_membership_data_term=0.0,
                    outmask_membership_data_term=(float(outmask_data[x, y, z]) if k == outmask_class else 0.0),
                )
                class_term = np.float32(update.class_membership_data_term)
                outmask_term = np.float32(update.outmask_membership_data_term)
                if coefficient > np.float32(0.0):
                    prior_update = update_lesion_toads_prior_energies(
                        atlas_priors,
                        similarity[k],
                        zero_similarities,
                        zero_similarities,
                        shapes[x, y, z],
                        class_index=k,
                        class_type=types[k],
                        prior_coefficient=float(coefficient),
                        class_membership_data_term=float(class_term),
                        lesion_membership_data_term=0.0,
                        black_hole_membership_data_term=0.0,
                        outmask_membership_data_term=float(outmask_term),
                        class_energy=0.0,
                        lesion_energy=0.0,
                        black_hole_energy=0.0,
                        has_lesion=False,
                        lesion_class_index=0,
                        lesion_class_type="OBJECT",
                    )
                    class_term = np.float32(prior_update.class_membership_data_term)
                    outmask_term = np.float32(prior_update.outmask_membership_data_term)
                class_membership = convert_lesion_toads_class_membership(
                    float(class_term),
                    float(shapes[x, y, z, k]),
                    relation_type=relations[k],
                    zero=float(zero_value),
                    infinity=float(infinity_value),
                )
                if k == outmask_class:
                    outmask_membership = convert_lesion_toads_class_membership(
                        float(outmask_term),
                        1.0,
                        relation_type="direct",
                        zero=float(zero_value),
                        infinity=float(infinity_value),
                    )
                    mixed = update_lesion_toads_outmask_membership(
                        float(convert_lesion_toads_class_membership(
                            float(class_term),
                            1.0,
                            relation_type="direct",
                            zero=float(zero_value),
                            infinity=float(infinity_value),
                        )),
                        float(outmask_membership),
                        float(outlier[x, y, z]),
                        relation_type=relations[k],
                        shape_prior=float(shapes[x, y, z, k]),
                    )
                    current[x, y, z, k] = mixed.class_membership
                    outlier[x, y, z] = mixed.outlier_membership
                else:
                    current[x, y, z, k] = class_membership
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if k == outmask_class:
                    outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)
    return LesionToadsOutmaskSweepResult(current, outlier)


def run_lesion_toads_lesion_outmask_smoothing_sweep(
    initial_class_memberships: np.ndarray,
    initial_lesion_membership: np.ndarray,
    initial_outlier_membership: np.ndarray,
    class_data_terms: np.ndarray,
    lesion_data_terms: np.ndarray,
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
) -> LesionToadsLesionOutmaskSweepResult:
    """Sweep OBJECT, LESION, and OUTMASK classes with shared auxiliaries.

    The function composes the selected lesion and outmask smoothing/prior
    terms in one source-ordered x/y/z pass. It updates the lesion and outlier
    maps in place during each voxel's class loop, then normalizes class and
    auxiliary memberships by the shared voxel denominator. Optional lesion
    distance maps apply after lesion inversion and before that normalization.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    outmask_data = np.asarray(outmask_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("LESION") != 1 or types.count("OUTMASK") != 1 or any(
        kind not in {"OBJECT", "LESION", "OUTMASK"} for kind in types
    ):
        raise ValueError("class_types must contain OBJECT classes and exactly one LESION and one OUTMASK class")
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    if any(array.shape != spatial_shape for array in (lesion, outlier, lesion_data, outmask_data)):
        raise ValueError("lesion/outlier memberships and auxiliary data maps must match class spatial shape")
    smoothing_priors_value = np.asarray(smoothing_priors, dtype=np.float32)
    if (
        smoothing_priors_value.shape != (classes,)
        or not np.isfinite(smoothing_priors_value).all()
        or np.any(smoothing_priors_value < 0)
    ):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not all(np.isfinite(array).all() for array in (current, lesion, outlier, data, lesion_data, outmask_data, shapes)):
        raise ValueError("membership, data-term, and shape arrays must be finite")
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient) or coefficient < np.float32(0.0):
        raise ValueError("prior_coefficient must be finite and nonnegative")
    atlas_priors = class_similarity = lesion_similarity = None
    if coefficient > np.float32(0.0):
        if any(value is None for value in (shape_priors, class_similarities, lesion_similarities)):
            raise ValueError("shape priors and class/lesion similarities are required for a positive prior_coefficient")
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        class_similarity = np.asarray(class_similarities, dtype=np.float32)
        lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or lesion_similarity.shape != (classes,):
            raise ValueError("shape_priors and lesion_similarities must have one value per class")
        if class_similarity.shape != (classes, classes):
            raise ValueError("class_similarities must have shape (classes, classes)")
        if not all(np.isfinite(array).all() for array in (atlas_priors, class_similarity, lesion_similarity)):
            raise ValueError("atlas prior and similarity inputs must be finite")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
    elif any(value is not None for value in (shape_priors, class_similarities, lesion_similarities)):
        if any(value is None for value in (shape_priors, class_similarities, lesion_similarities)):
            raise ValueError("shape_priors and class/lesion similarities must be provided together")
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        class_similarity = np.asarray(class_similarities, dtype=np.float32)
        lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or lesion_similarity.shape != (classes,):
            raise ValueError("shape_priors and lesion_similarities must have one value per class")
        if class_similarity.shape != (classes, classes):
            raise ValueError("class_similarities must have shape (classes, classes)")
        if not all(np.isfinite(array).all() for array in (atlas_priors, class_similarity, lesion_similarity)):
            raise ValueError("atlas prior and similarity inputs must be finite")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
    maxima = (
        ("max_ventricle_distance", max_ventricle_distance),
        ("max_gray_matter_distance", max_gray_matter_distance),
        ("max_brainstem_distance", max_brainstem_distance),
    )
    for name, value in maxima:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not -128 <= int(value) <= 127:
            raise ValueError(f"{name} must fit a signed byte")
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")
    spread_value = np.float32(spread)
    if not np.isfinite(spread_value):
        raise ValueError("spread must be finite")
    distance_inputs = (
        ("ventricle_distances", ventricle_distances, max_ventricle_distance),
        ("gray_matter_distances", gray_matter_distances, max_gray_matter_distance),
        ("brainstem_distances", brainstem_distances, max_brainstem_distance),
    )
    distance_maps: dict[str, np.ndarray] = {}
    for name, source, maximum in distance_inputs:
        if source is None:
            if maximum != 0:
                raise ValueError(f"{name} is required when its maximum distance is nonzero")
            distance_maps[name] = np.full(spatial_shape, -1, dtype=np.int8)
            continue
        raw = np.asarray(source)
        if not np.issubdtype(raw.dtype, np.integer) or raw.shape != spatial_shape:
            raise ValueError(f"{name} must be an integer map matching the class spatial shape")
        if np.any(raw < -1) or np.any(raw > 127):
            raise ValueError(f"{name} values must fit the source -1..127 distance range")
        distance_maps[name] = np.asarray(raw, dtype=np.int8)
    use_distances = any(value != 0 for _, value in maxima)
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    zero_map = np.zeros(spatial_shape, dtype=np.float32)
    zero_similarity = np.zeros(classes, dtype=np.float32)
    lesion_class = types.index("LESION")
    outmask_class = types.index("OUTMASK")
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, zero_map):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            if shapes[x, y, z, k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
                if k == lesion_class:
                    lesion[x, y, z] = np.float32(0.0)
            else:
                update = update_lesion_toads_smoothing_energies(
                    neighborhood.class_memberships,
                    neighborhood.lesion_memberships,
                    neighborhood.black_hole_memberships,
                    class_index=k,
                    class_type=types[k],
                    smoothing_prior=float(smoothing_priors_value[k]),
                    smoothing=float(smoothing_value),
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    class_membership_data_term=float(data[x, y, z, k]),
                    lesion_membership_data_term=(float(lesion_data[x, y, z]) if k == lesion_class else 0.0),
                    black_hole_membership_data_term=0.0,
                    outmask_membership_data_term=(float(outmask_data[x, y, z]) if k == outmask_class else 0.0),
                )
                class_term = np.float32(update.class_membership_data_term)
                lesion_term = np.float32(update.lesion_membership_data_term)
                outmask_term = np.float32(update.outmask_membership_data_term)
                if coefficient > np.float32(0.0):
                    prior_update = update_lesion_toads_prior_energies(
                        atlas_priors,
                        class_similarity[k],
                        lesion_similarity,
                        zero_similarity,
                        shapes[x, y, z],
                        class_index=k,
                        class_type=types[k],
                        prior_coefficient=float(coefficient),
                        class_membership_data_term=float(class_term),
                        lesion_membership_data_term=float(lesion_term),
                        black_hole_membership_data_term=0.0,
                        outmask_membership_data_term=float(outmask_term),
                        class_energy=0.0,
                        lesion_energy=0.0,
                        black_hole_energy=0.0,
                        has_lesion=True,
                        lesion_class_index=lesion_class,
                        lesion_class_type="LESION",
                    )
                    class_term = np.float32(prior_update.class_membership_data_term)
                    lesion_term = np.float32(prior_update.lesion_membership_data_term)
                    outmask_term = np.float32(prior_update.outmask_membership_data_term)
                class_membership = convert_lesion_toads_class_membership(
                    float(class_term), float(shapes[x, y, z, k]), relation_type=relations[k],
                    zero=float(zero_value), infinity=float(infinity_value),
                )
                if k == lesion_class:
                    lesion_membership = convert_lesion_toads_class_membership(
                        float(lesion_term), float(shapes[x, y, z, k]), relation_type=relations[k],
                        zero=float(zero_value), infinity=float(infinity_value),
                    )
                    lesion[x, y, z] = lesion_membership
                    if use_distances:
                        shape = np.float32(shapes[x, y, z, k])
                        lesion_relation = (
                            np.float32(1.0)
                            if relations[k] == "direct"
                            else np.float32(1.0)
                            if np.float32(1.0) < np.float32(np.float32(2.0) * shape)
                            else np.float32(np.float32(2.0) * shape)
                        )
                        modulated = modulate_lesion_toads_membership_by_distances(
                            float(lesion[x, y, z]),
                            float(lesion_relation),
                            ventricle_distance=int(distance_maps["ventricle_distances"][x, y, z]),
                            gray_matter_distance=int(distance_maps["gray_matter_distances"][x, y, z]),
                            brainstem_distance=int(distance_maps["brainstem_distances"][x, y, z]),
                            max_ventricle_distance=int(max_ventricle_distance),
                            max_gray_matter_distance=int(max_gray_matter_distance),
                            max_brainstem_distance=int(max_brainstem_distance),
                            gm_flag=gm_flag,
                            spread=float(spread_value),
                        )
                        lesion[x, y, z] = modulated.membership
                if k == outmask_class:
                    outmask_membership = convert_lesion_toads_class_membership(
                        float(outmask_term), 1.0, relation_type="direct",
                        zero=float(zero_value), infinity=float(infinity_value),
                    )
                    mixed = update_lesion_toads_outmask_membership(
                        float(convert_lesion_toads_class_membership(
                            float(class_term), 1.0, relation_type="direct",
                            zero=float(zero_value), infinity=float(infinity_value),
                        )),
                        float(outmask_membership), float(outlier[x, y, z]),
                        relation_type=relations[k], shape_prior=float(shapes[x, y, z, k]),
                    )
                    current[x, y, z, k] = mixed.class_membership
                    outlier[x, y, z] = mixed.outlier_membership
                else:
                    current[x, y, z, k] = class_membership
            if k == lesion_class:
                denominator = np.float32(denominator + lesion[x, y, z])
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if k == lesion_class:
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                elif k == outmask_class:
                    outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)
    return LesionToadsLesionOutmaskSweepResult(current, lesion, outlier)


def run_lesion_toads_blackhole_outmask_smoothing_sweep(
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
    """Sweep OBJECT, optional LESION, BLACKHOLE, and OUTMASK classes.

    This composes the selected class, lesion, black-hole, and OUTMASK energy
    updates in source order, applies BLACKHOLE distance rules after inversion,
    then normalizes all memberships by the shared voxel denominator.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    black_hole = np.asarray(initial_black_hole_membership, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    black_hole_data = np.asarray(black_hole_data_terms, dtype=np.float32)
    outmask_data = np.asarray(outmask_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    if current.ndim != 4 or any(size < 3 for size in current.shape[:3]) or current.shape[3] == 0:
        raise ValueError("initial_class_memberships must have shape (x, y, z, classes), with spatial sizes >= 3")
    spatial_shape = current.shape[:3]
    classes = current.shape[3]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("BLACKHOLE") != 1 or types.count("OUTMASK") != 1 or types.count("LESION") > 1 or any(
        kind not in {"OBJECT", "LESION", "BLACKHOLE", "OUTMASK"} for kind in types
    ):
        raise ValueError("class_types must contain OBJECT classes, exactly one BLACKHOLE and one OUTMASK, and at most one LESION")
    if len(relations) != classes or any(value not in {"direct", "prior"} for value in relations):
        raise ValueError("relation_types must contain one 'direct' or 'prior' value per class")
    if data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match initial_class_memberships shape")
    if any(array.shape != spatial_shape for array in (lesion, black_hole, outlier, lesion_data, black_hole_data, outmask_data)):
        raise ValueError("auxiliary membership and data maps must match the class spatial shape")
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    if smooth_prior.shape != (classes,) or not np.isfinite(smooth_prior).all() or np.any(smooth_prior < 0):
        raise ValueError("smoothing_priors must be finite nonnegative values, one per class")
    if not all(np.isfinite(array).all() for array in (current, lesion, black_hole, outlier, data, lesion_data, black_hole_data, outmask_data, shapes)):
        raise ValueError("membership, data-term, and shape arrays must be finite")
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient) or coefficient < np.float32(0.0):
        raise ValueError("prior_coefficient must be finite and nonnegative")
    atlas_priors = class_similarity = lesion_similarity = black_hole_similarity = None
    prior_inputs = (shape_priors, class_similarities, lesion_similarities, black_hole_similarities)
    if coefficient > np.float32(0.0) and any(value is None for value in prior_inputs):
        raise ValueError("shape priors and all class/auxiliary similarities are required for a positive prior_coefficient")
    if any(value is not None for value in prior_inputs):
        if any(value is None for value in prior_inputs):
            raise ValueError("shape priors and all class/auxiliary similarities must be provided together")
        atlas_priors = np.asarray(shape_priors, dtype=np.float32)
        class_similarity = np.asarray(class_similarities, dtype=np.float32)
        lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
        black_hole_similarity = np.asarray(black_hole_similarities, dtype=np.float32)
        if atlas_priors.shape != (classes,) or lesion_similarity.shape != (classes,) or black_hole_similarity.shape != (classes,):
            raise ValueError("shape priors and auxiliary similarities must have one value per class")
        if class_similarity.shape != (classes, classes):
            raise ValueError("class_similarities must have shape (classes, classes)")
        if not all(np.isfinite(array).all() for array in (atlas_priors, class_similarity, lesion_similarity, black_hole_similarity)):
            raise ValueError("atlas prior and similarity inputs must be finite")
        if np.any(atlas_priors < np.float32(0.0)):
            raise ValueError("shape_priors must be nonnegative")
    smoothing_value = np.float32(smoothing)
    if not np.isfinite(smoothing_value) or smoothing_value <= 0:
        raise ValueError("smoothing must be finite and positive")
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(zero_value) or zero_value < 0:
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("infinity must be finite and nonnegative")

    maxima = (
        ("max_ventricle_distance", max_ventricle_distance),
        ("max_gray_matter_distance", max_gray_matter_distance),
        ("max_brainstem_distance", max_brainstem_distance),
    )
    for name, value in maxima:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not -128 <= int(value) <= 127:
            raise ValueError(f"{name} must fit a signed byte")
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")
    spread_value = np.float32(spread)
    if not np.isfinite(spread_value):
        raise ValueError("spread must be finite")
    distance_inputs = (
        ("ventricle_distances", ventricle_distances, max_ventricle_distance),
        ("gray_matter_distances", gray_matter_distances, max_gray_matter_distance),
        ("brainstem_distances", brainstem_distances, max_brainstem_distance),
    )
    distance_maps: dict[str, np.ndarray] = {}
    for name, source, maximum in distance_inputs:
        if source is None:
            if maximum != 0:
                raise ValueError(f"{name} is required when its maximum distance is nonzero")
            distance_maps[name] = np.full(spatial_shape, -1, dtype=np.int8)
            continue
        raw = np.asarray(source)
        if not np.issubdtype(raw.dtype, np.integer) or raw.shape != spatial_shape:
            raise ValueError(f"{name} must be an integer map matching the class spatial shape")
        if np.any(raw < -1) or np.any(raw > 127):
            raise ValueError(f"{name} values must fit the source -1..127 distance range")
        distance_maps[name] = np.asarray(raw, dtype=np.int8)
    use_distances = any(value != 0 for _, value in maxima)

    black_hole_class = types.index("BLACKHOLE")
    outmask_class = types.index("OUTMASK")
    lesion_class = max((k for k, kind in enumerate(types) if kind in {"LESION", "BLACKHOLE"}), default=black_hole_class)
    lesion_class_type = types[lesion_class]
    zero_similarity = np.zeros(classes, dtype=np.float32)
    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k in range(classes):
            if shapes[x, y, z, k] <= np.float32(0.0):
                current[x, y, z, k] = np.float32(0.0)
                if types[k] in {"LESION", "BLACKHOLE"}:
                    lesion[x, y, z] = np.float32(0.0)
                if k == black_hole_class:
                    black_hole[x, y, z] = np.float32(0.0)
            else:
                update = update_lesion_toads_smoothing_energies(
                    neighborhood.class_memberships,
                    neighborhood.lesion_memberships,
                    neighborhood.black_hole_memberships,
                    class_index=k,
                    class_type=types[k],
                    smoothing_prior=float(smooth_prior[k]),
                    smoothing=float(smoothing_value),
                    class_energy=0.0,
                    lesion_energy=0.0,
                    black_hole_energy=0.0,
                    class_membership_data_term=float(data[x, y, z, k]),
                    lesion_membership_data_term=(float(lesion_data[x, y, z]) if types[k] in {"LESION", "BLACKHOLE"} else 0.0),
                    black_hole_membership_data_term=(float(black_hole_data[x, y, z]) if k == black_hole_class else 0.0),
                    outmask_membership_data_term=(float(outmask_data[x, y, z]) if k == outmask_class else 0.0),
                )
                class_term = np.float32(update.class_membership_data_term)
                lesion_term = np.float32(update.lesion_membership_data_term)
                black_hole_term = np.float32(update.black_hole_membership_data_term)
                outmask_term = np.float32(update.outmask_membership_data_term)
                if coefficient > np.float32(0.0):
                    prior_update = update_lesion_toads_prior_energies(
                        atlas_priors,
                        class_similarity[k],
                        lesion_similarity,
                        black_hole_similarity,
                        shapes[x, y, z],
                        class_index=k,
                        class_type=types[k],
                        prior_coefficient=float(coefficient),
                        class_membership_data_term=float(class_term),
                        lesion_membership_data_term=float(lesion_term),
                        black_hole_membership_data_term=float(black_hole_term),
                        outmask_membership_data_term=float(outmask_term),
                        class_energy=0.0,
                        lesion_energy=0.0,
                        black_hole_energy=0.0,
                        has_lesion=True,
                        lesion_class_index=lesion_class,
                        lesion_class_type=lesion_class_type,
                    )
                    class_term = np.float32(prior_update.class_membership_data_term)
                    lesion_term = np.float32(prior_update.lesion_membership_data_term)
                    black_hole_term = np.float32(prior_update.black_hole_membership_data_term)
                    outmask_term = np.float32(prior_update.outmask_membership_data_term)
                current[x, y, z, k] = convert_lesion_toads_class_membership(
                    float(class_term), float(shapes[x, y, z, k]), relation_type=relations[k],
                    zero=float(zero_value), infinity=float(infinity_value),
                )
                if types[k] in {"LESION", "BLACKHOLE"}:
                    lesion[x, y, z] = convert_lesion_toads_class_membership(
                        float(lesion_term), float(shapes[x, y, z, k]), relation_type=relations[k],
                        zero=float(zero_value), infinity=float(infinity_value),
                    )
                if k == black_hole_class:
                    black_hole[x, y, z] = convert_lesion_toads_class_membership(
                        float(black_hole_term), float(shapes[x, y, z, k]), relation_type=relations[k],
                        zero=float(zero_value), infinity=float(infinity_value),
                    )
                    if use_distances:
                        shape = np.float32(shapes[x, y, z, k])
                        lesion_relation = (
                            np.float32(1.0) if relations[k] == "direct"
                            else np.float32(1.0) if np.float32(1.0) < np.float32(np.float32(2.0) * shape)
                            else np.float32(np.float32(2.0) * shape)
                        )
                        modulated = modulate_lesion_toads_blackhole_memberships_by_distances(
                            float(lesion[x, y, z]), float(black_hole[x, y, z]), float(lesion_relation),
                            ventricle_distance=int(distance_maps["ventricle_distances"][x, y, z]),
                            gray_matter_distance=int(distance_maps["gray_matter_distances"][x, y, z]),
                            brainstem_distance=int(distance_maps["brainstem_distances"][x, y, z]),
                            max_ventricle_distance=int(max_ventricle_distance),
                            max_gray_matter_distance=int(max_gray_matter_distance),
                            max_brainstem_distance=int(max_brainstem_distance),
                            gm_flag=gm_flag, spread=float(spread_value),
                        )
                        lesion[x, y, z] = modulated.lesion_membership
                        black_hole[x, y, z] = modulated.black_hole_membership
                elif types[k] == "LESION" and use_distances:
                    shape = np.float32(shapes[x, y, z, k])
                    lesion_relation = (
                        np.float32(1.0) if relations[k] == "direct"
                        else np.float32(1.0) if np.float32(1.0) < np.float32(np.float32(2.0) * shape)
                        else np.float32(np.float32(2.0) * shape)
                    )
                    modulated = modulate_lesion_toads_membership_by_distances(
                        float(lesion[x, y, z]), float(lesion_relation),
                        ventricle_distance=int(distance_maps["ventricle_distances"][x, y, z]),
                        gray_matter_distance=int(distance_maps["gray_matter_distances"][x, y, z]),
                        brainstem_distance=int(distance_maps["brainstem_distances"][x, y, z]),
                        max_ventricle_distance=int(max_ventricle_distance),
                        max_gray_matter_distance=int(max_gray_matter_distance),
                        max_brainstem_distance=int(max_brainstem_distance),
                        gm_flag=gm_flag, spread=float(spread_value),
                    )
                    lesion[x, y, z] = modulated.membership
                if k == outmask_class:
                    outmask_membership = convert_lesion_toads_class_membership(
                        float(outmask_term), 1.0, relation_type="direct",
                        zero=float(zero_value), infinity=float(infinity_value),
                    )
                    mixed = update_lesion_toads_outmask_membership(
                        float(convert_lesion_toads_class_membership(
                            float(class_term), 1.0, relation_type="direct",
                            zero=float(zero_value), infinity=float(infinity_value),
                        )),
                        float(outmask_membership), float(outlier[x, y, z]),
                        relation_type=relations[k], shape_prior=float(shapes[x, y, z, k]),
                    )
                    current[x, y, z, k] = mixed.class_membership
                    outlier[x, y, z] = mixed.outlier_membership
            if types[k] == "LESION":
                denominator = np.float32(denominator + lesion[x, y, z])
            elif k == black_hole_class:
                auxiliary_sum = np.float32(lesion[x, y, z] + black_hole[x, y, z])
                denominator = np.float32(denominator + auxiliary_sum)
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > np.float32(0.0):
            for k in range(classes):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if types[k] in {"LESION", "BLACKHOLE"}:
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                if k == black_hole_class:
                    black_hole[x, y, z] = np.float32(black_hole[x, y, z] / denominator)
                elif k == outmask_class:
                    outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)
    return LesionToadsBlackHoleOutmaskSweepResult(current, lesion, black_hole, outlier)


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


def update_lesion_toads_smoothing_energies(
    class_neighbor_memberships: np.ndarray,
    lesion_neighbor_memberships: np.ndarray,
    black_hole_neighbor_memberships: np.ndarray,
    *,
    class_index: int,
    class_type: str,
    smoothing_prior: float,
    smoothing: float,
    class_energy: float,
    lesion_energy: float,
    black_hole_energy: float,
    class_membership_data_term: float,
    lesion_membership_data_term: float,
    black_hole_membership_data_term: float,
    outmask_membership_data_term: float,
) -> LesionToadsSmoothingEnergyUpdate:
    """Apply one source-ordered six-face smoothing update at a voxel.

    Face order for the input is X+, X-, Y+, Y-, Z+, Z-. The selected source
    accumulates all other-class neighbor squares in that face order and class
    order, then separately sums the current class's self-neighbor term in
    X-, X+, Y-, Y+, Z-, Z+ order for lesion-related types.
    """
    memberships = np.asarray(class_neighbor_memberships, dtype=np.float32)
    lesion_faces = np.asarray(lesion_neighbor_memberships, dtype=np.float32)
    black_hole_faces = np.asarray(black_hole_neighbor_memberships, dtype=np.float32)
    if memberships.ndim != 2 or memberships.shape[0] != 6 or memberships.shape[1] == 0:
        raise ValueError("class_neighbor_memberships must have shape (6, classes)")
    if lesion_faces.shape != (6,) or black_hole_faces.shape != (6,):
        raise ValueError("auxiliary neighbor memberships must each have six face values")
    if isinstance(class_index, (bool, np.bool_)) or not isinstance(class_index, (int, np.integer)):
        raise ValueError("class_index must be an integer")
    if not 0 <= int(class_index) < memberships.shape[1]:
        raise ValueError("class_index must identify one input membership channel")
    supported = {"OBJECT", "LESION", "OPTIMIZED_LESION", "BLACKHOLE", "OUTMASK", "OUTLIER", "MASK"}
    if class_type not in supported:
        raise ValueError("class_type must be OBJECT, LESION, OPTIMIZED_LESION, BLACKHOLE, or OUTMASK")

    k = int(class_index)
    other = np.float32(0.0)
    lesion_sum = np.float32(0.0)
    black_hole_sum = np.float32(0.0)
    # The Java source visits X+, X-, Y+, Y-, Z+, Z- and each class in order.
    for face in range(6):
        for class_no in range(memberships.shape[1]):
            if class_no != k:
                square = np.float32(memberships[face, class_no] * memberships[face, class_no])
                other = np.float32(other + square)
        lesion_square = np.float32(lesion_faces[face] * lesion_faces[face])
        lesion_sum = np.float32(lesion_sum + lesion_square)
        black_hole_square = np.float32(black_hole_faces[face] * black_hole_faces[face])
        black_hole_sum = np.float32(black_hole_sum + black_hole_square)

    current_sum = np.float32(0.0)
    if class_type in {"LESION", "OPTIMIZED_LESION", "BLACKHOLE"}:
        for face in (1, 0, 3, 2, 5, 4):
            square = np.float32(memberships[face, k] * memberships[face, k])
            current_sum = np.float32(current_sum + square)

    p = np.float32(smoothing_prior)
    smooth = np.float32(smoothing)
    cls_membership_term = np.float32(class_membership_data_term)
    lesion_membership_term = np.float32(lesion_membership_data_term)
    black_hole_membership_term = np.float32(black_hole_membership_data_term)
    outmask_membership_term = np.float32(outmask_membership_data_term)
    cls_energy = np.float32(class_energy)
    lesion_e = np.float32(lesion_energy)
    black_hole_e = np.float32(black_hole_energy)

    all_neighbor_sum = np.float32(other + lesion_sum)
    all_neighbor_sum = np.float32(all_neighbor_sum + black_hole_sum)
    class_full = np.float32(np.float32(p * smooth) * all_neighbor_sum)
    class_half = np.float32(np.float32(np.float32(np.float32(0.5) * p) * smooth) * all_neighbor_sum)
    cls_membership_term = np.float32(cls_membership_term + class_full)
    cls_energy = np.float32(cls_energy + class_half)

    if class_type in {"LESION", "OPTIMIZED_LESION", "BLACKHOLE"}:
        lesion_sum_for_update = np.float32(other + current_sum)
        lesion_sum_for_update = np.float32(lesion_sum_for_update + black_hole_sum)
        lesion_full = np.float32(np.float32(p * smooth) * lesion_sum_for_update)
        lesion_half = np.float32(np.float32(np.float32(np.float32(0.5) * p) * smooth) * lesion_sum_for_update)
        lesion_membership_term = np.float32(lesion_membership_term + lesion_full)
        lesion_e = np.float32(lesion_e + lesion_half)

    if class_type == "BLACKHOLE":
        black_hole_sum_for_update = np.float32(other + current_sum)
        black_hole_sum_for_update = np.float32(black_hole_sum_for_update + lesion_sum)
        black_hole_full = np.float32(np.float32(p * smooth) * black_hole_sum_for_update)
        black_hole_half = np.float32(
            np.float32(np.float32(np.float32(np.float32(0.5) * p) * smooth) * black_hole_sum_for_update)
        )
        black_hole_membership_term = np.float32(black_hole_membership_term + black_hole_full)
        black_hole_e = np.float32(black_hole_e + black_hole_half)

    if class_type == "OUTMASK":
        outmask_sum = np.float32(other + lesion_sum)
        outmask_full = np.float32(np.float32(p * smooth) * outmask_sum)
        outmask_membership_term = np.float32(outmask_membership_term + outmask_full)

    return LesionToadsSmoothingEnergyUpdate(
        class_membership_data_term=cls_membership_term,
        lesion_membership_data_term=lesion_membership_term,
        black_hole_membership_data_term=black_hole_membership_term,
        outmask_membership_data_term=outmask_membership_term,
        class_energy=cls_energy,
        lesion_energy=lesion_e,
        black_hole_energy=black_hole_e,
        other_class_neighbor_squares=other,
        lesion_neighbor_squares=lesion_sum,
        black_hole_neighbor_squares=black_hole_sum,
        current_class_neighbor_squares=current_sum,
    )


def update_lesion_toads_prior_energies(
    shape_priors: np.ndarray,
    class_similarities: np.ndarray,
    lesion_similarities: np.ndarray,
    black_hole_similarities: np.ndarray,
    shape_values: np.ndarray,
    *,
    class_index: int,
    class_type: str,
    prior_coefficient: float,
    class_membership_data_term: float,
    lesion_membership_data_term: float,
    black_hole_membership_data_term: float,
    outmask_membership_data_term: float,
    class_energy: float,
    lesion_energy: float,
    black_hole_energy: float,
    has_lesion: bool,
    lesion_class_index: int,
    lesion_class_type: str,
) -> LesionToadsPriorEnergyUpdate:
    """Apply one selected prior-coefficient update at one voxel.

    Similarity vectors are the selected row similarity[k][m],
    lesionSimilarity[m], and blackHoleSimilarity[m]. Products and sums
    follow the Java loop order and float32 left-associative arithmetic.
    """
    priors = np.asarray(shape_priors, dtype=np.float32)
    similarities = np.asarray(class_similarities, dtype=np.float32)
    lesion_sim = np.asarray(lesion_similarities, dtype=np.float32)
    black_hole_sim = np.asarray(black_hole_similarities, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    arrays = (priors, similarities, lesion_sim, black_hole_sim, shapes)
    if any(array.ndim != 1 for array in arrays):
        raise ValueError("prior and similarity inputs must be vectors")
    if not priors.size or any(array.size != priors.size for array in arrays[1:]):
        raise ValueError("prior and similarity vectors must have equal nonzero length")
    if isinstance(class_index, (bool, np.bool_)) or not isinstance(class_index, (int, np.integer)):
        raise ValueError("class_index must be an integer")
    if not 0 <= int(class_index) < priors.size:
        raise ValueError("class_index must identify an input class")
    if not isinstance(has_lesion, (bool, np.bool_)):
        raise ValueError("has_lesion must be boolean")
    if has_lesion and (
        isinstance(lesion_class_index, (bool, np.bool_))
        or not isinstance(lesion_class_index, (int, np.integer))
        or not 0 <= int(lesion_class_index) < priors.size
    ):
        raise ValueError("lesion_class_index must identify a class when has_lesion is true")
    if class_type not in {"OBJECT", "LESION", "OPTIMIZED_LESION", "BLACKHOLE", "OUTMASK", "OUTLIER", "MASK"}:
        raise ValueError("unsupported LesionTOADS class_type")

    k = int(class_index)
    class_prior = np.float32(0.0)
    lesion_prior = np.float32(0.0)
    black_hole_prior = np.float32(0.0)
    for m in range(priors.size):
        shape = shapes[m]
        if m != k:
            product = np.float32(priors[m] * similarities[m])
            product = np.float32(product * shape)
            product = np.float32(product * shape)
            class_prior = np.float32(class_prior + product)
        if class_type in {"LESION", "OPTIMIZED_LESION", "BLACKHOLE"}:
            product = np.float32(priors[m] * lesion_sim[m])
            product = np.float32(product * shape)
            product = np.float32(product * shape)
            lesion_prior = np.float32(lesion_prior + product)
        if class_type == "BLACKHOLE":
            product = np.float32(priors[m] * black_hole_sim[m])
            product = np.float32(product * shape)
            product = np.float32(product * shape)
            black_hole_prior = np.float32(black_hole_prior + product)

    if has_lesion:
        lesion_index = int(lesion_class_index)
        shape = shapes[lesion_index]
        product = np.float32(priors[lesion_index] * lesion_sim[k])
        product = np.float32(product * shape)
        product = np.float32(product * shape)
        class_prior = np.float32(class_prior + product)
        if lesion_class_type == "BLACKHOLE":
            product = np.float32(priors[lesion_index] * black_hole_sim[k])
            product = np.float32(product * shape)
            product = np.float32(product * shape)
            class_prior = np.float32(class_prior + product)

    coefficient = np.float32(prior_coefficient)
    class_term = np.float32(class_membership_data_term)
    lesion_term = np.float32(lesion_membership_data_term)
    black_hole_term = np.float32(black_hole_membership_data_term)
    outmask_term = np.float32(outmask_membership_data_term)
    class_term = np.float32(class_term + np.float32(coefficient * class_prior))
    if class_type == "OUTMASK":
        outmask_term = np.float32(outmask_term + np.float32(coefficient * class_prior))
    class_energy_value = np.float32(class_energy)
    class_half = np.float32(np.float32(np.float32(0.5) * coefficient) * class_prior)
    class_energy_value = np.float32(class_energy_value + class_half)

    lesion_energy_value = np.float32(lesion_energy)
    if class_type in {"LESION", "OPTIMIZED_LESION", "BLACKHOLE"}:
        lesion_term = np.float32(lesion_term + np.float32(coefficient * lesion_prior))
        half = np.float32(np.float32(np.float32(0.5) * coefficient) * lesion_prior)
        lesion_energy_value = np.float32(lesion_energy_value + half)

    black_hole_energy_value = np.float32(black_hole_energy)
    if class_type == "BLACKHOLE":
        black_hole_term = np.float32(black_hole_term + np.float32(coefficient * black_hole_prior))
        half = np.float32(np.float32(np.float32(0.5) * coefficient) * black_hole_prior)
        black_hole_energy_value = np.float32(black_hole_energy_value + half)

    return LesionToadsPriorEnergyUpdate(
        class_prior=class_prior,
        lesion_prior=lesion_prior,
        black_hole_prior=black_hole_prior,
        class_membership_data_term=class_term,
        lesion_membership_data_term=lesion_term,
        black_hole_membership_data_term=black_hole_term,
        outmask_membership_data_term=outmask_term,
        class_energy=class_energy_value,
        lesion_energy=lesion_energy_value,
        black_hole_energy=black_hole_energy_value,
    )


def update_lesion_toads_outmask_membership(
    class_membership: float,
    outmask_membership: float,
    outlier_membership: float,
    *,
    relation_type: str,
    shape_prior: float,
) -> LesionToadsOutmaskMembership:
    """Blend the selected OUTMASK and class memberships by outlier weight.

    This is the `OUTMASK` assignment after the class and outmask terms have
    already been inverted. Subsequent shared denominator normalization is a
    separate source step.
    """
    class_value = np.float32(class_membership)
    outmask_value = np.float32(outmask_membership)
    outlier = np.float32(outlier_membership)
    shape = np.float32(shape_prior)
    if not all(np.isfinite(value) for value in (class_value, outmask_value, outlier, shape)):
        raise ValueError("membership and shape inputs must be finite")
    if relation_type not in {"direct", "prior"}:
        raise ValueError("relation_type must be 'direct' or 'prior'")

    class_fraction = np.float32(np.float32(1.0) - outlier)
    class_component = np.float32(class_fraction * class_value)
    outlier_component = np.float32(outlier * outmask_value)
    mixed = np.float32(class_component + outlier_component)
    if relation_type == "prior":
        doubled_shape = np.float32(np.float32(2.0) * shape)
        shape_factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
        mixed = np.float32(shape_factor * mixed)
        outlier = np.float32(shape_factor * outlier)
    outlier = np.float32(outlier * outmask_value)
    return LesionToadsOutmaskMembership(
        class_membership=mixed,
        outlier_membership=outlier,
    )
