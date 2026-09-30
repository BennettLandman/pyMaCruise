"""Fixed-input exact-membership dispatch for selected atlas topology types.

This composes the bundled TOADS2012 13-structure descriptor with the active
``MultipleToadSegmentation.computeExactMemberships`` sweep. Atlas statistics,
transformed subject-grid priors, and relation maps are explicit inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_object_intensity_energy,
    toads_optimized_tissue_data_term,
    toads_shape_relation_factors,
)


_TOPOLOGY_TYPE = {
    "obj": "OBJECT",
    "mask": "MASK",
    "out": "OUTLIER",
    "lesion": "LESION",
    "outmask": "OUTMASK",
    "optimized_csf": "OPTIMIZED_CSF",
    "optimized_gm": "OPTIMIZED_GM",
    "optimized_wm": "OPTIMIZED_WM",
    "rayleigh": "RAYLEIGH",
}
_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}


@dataclass(frozen=True)
class ToadsAtlasTypeMembershipPass:
    """Class energies, memberships, denominator, and objective for one pass."""

    class_types: tuple[str, ...]
    class_data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def map_toads_topologies(topologies: Sequence[str]) -> tuple[str, ...]:
    """Map exact selected atlas topology tokens to segmentation class types.

    The Java constructor uses case-sensitive string equality and falls back to
    OBJECT for an unrecognized token. Preserve both behaviors.
    """
    return tuple(_TOPOLOGY_TYPE.get(str(value), "OBJECT") for value in topologies)


def read_toads_atlas_class_types(path: str | Path) -> tuple[str, ...]:
    """Read class topology tokens from a selected Structure Atlas text file.

    This mirrors the descriptor's ordered ``Structures`` block; it does not
    load image payloads or attempt to interpret other atlas sections.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "Structure Atlas File (edit at your own risks)":
        raise ValueError("not a selected Structure Atlas file")
    header_index = next((i for i, line in enumerate(lines) if line.startswith("Structures:")), None)
    if header_index is None:
        raise ValueError("Structure Atlas file has no Structures block")
    try:
        count = int(lines[header_index].split(":", 1)[1].strip())
    except (IndexError, ValueError) as error:
        raise ValueError("invalid Structures count") from error
    if count < 1 or header_index + count >= len(lines):
        raise ValueError("truncated Structures block")
    topology_tokens: list[str] = []
    for line in lines[header_index + 1 : header_index + 1 + count]:
        columns = line.split("\t")
        if len(columns) < 3 or not columns[0]:
            raise ValueError("malformed Structure Atlas class row")
        topology_tokens.append(columns[2])
    return map_toads_topologies(topology_tokens)


def toads_atlas_type_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    class_types: Sequence[str],
    centroids: ArrayLike,
    exponent_by_channel_class: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    transformed_shapes: ArrayLike,
    closest_labels: ArrayLike,
    smoothing_priors: ArrayLike,
    centroid_similarity: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsAtlasTypeMembershipPass:
    """Compose exact membership updates for OBJECT/MASK/optimized classes.

    Inputs are fixed model state on the subject's true-RAS voxel grid. The
    selected atlas types are class ordered. This bounded adapter supports the
    actual 13-class TOADS2012 mix (MASK, optimized CSF, and OBJECT); unsupported
    active class families are rejected explicitly.
    """
    types = tuple(str(value).upper() for value in class_types)
    if not types or types.count("MASK") != 1 or types.count("OPTIMIZED_CSF") != 2:
        raise ValueError("class_types must describe the selected 13-class TOADS2012 topology order")
    if any(kind not in {"OBJECT", "MASK", "OPTIMIZED_CSF"} for kind in types):
        raise ValueError("class_types contains unsupported families; this adapter supports OBJECT, MASK, and OPTIMIZED_CSF")
    if len(types) != 13 or types[0] != "MASK" or types[1:3] != ("OPTIMIZED_CSF",) * 2:
        raise ValueError("class_types must match the selected 13-class TOADS2012 order")

    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    var_weights = np.asarray(class_variance_weights, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    labels = np.asarray(closest_labels)
    if images.ndim != 4 or images.shape[-1] < 1 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must have shape (x>=3,y>=3,z>=3,channels)")
    spatial = images.shape[:3]
    channels, classes = images.shape[-1], len(types)
    if initial.shape != spatial + (classes,) or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must match the 13-class subject grid")
    if centers.shape != (channels, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite with shape (channels,13)")
    if exponents.shape != centers.shape or not np.isfinite(exponents).all():
        raise ValueError("exponent_by_channel_class must be finite with shape (channels,13)")
    if var_weights.shape != centers.shape or not np.isfinite(var_weights).all():
        raise ValueError("class_variance_weights must be finite with shape (channels,13)")
    if shapes.shape != initial.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match subject memberships")
    if labels.shape != spatial + (4,) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot subject-grid map")

    modality = np.asarray(modality_weights, dtype=np.float32)
    image_var = np.asarray(image_variances, dtype=np.float32)
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    similarity = np.asarray(centroid_similarity, dtype=np.float32)
    shape_prior = np.asarray(atlas_shape_priors, dtype=np.float32)
    if modality.shape != (channels,) or image_var.shape != (channels,):
        raise ValueError("modality_weights and image_variances must match channel count")
    if smooth_prior.shape != (classes,) or not np.isfinite(smooth_prior).all():
        raise ValueError("smoothing_priors must contain 13 finite values")
    if similarity.shape != (classes, classes) or not np.isfinite(similarity).all():
        raise ValueError("centroid_similarity must be finite and 13 by 13")
    if shape_prior.shape != (classes,) or not np.isfinite(shape_prior).all():
        raise ValueError("atlas_shape_priors must contain 13 finite values")
    if not np.isfinite(images).all() or not np.isfinite(modality).all() or not np.isfinite(image_var).all():
        raise ValueError("image and intensity weights must be finite")

    for kind in _OPTIMIZED_TYPES:
        columns = [index for index, value in enumerate(types) if value == kind]
        if len(columns) > 1 and not np.array_equal(
            exponents[:, columns],
            np.broadcast_to(exponents[:, columns[0], None], (channels, len(columns))),
        ):
            raise ValueError(f"all {kind} classes must share category/channel exponents")

    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    bias = None if psf_bias is None else np.asarray(psf_bias, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")
    if bias is not None and (bias.shape != images.shape or not np.isfinite(bias).all()):
        raise ValueError("psf_bias must be finite and match image_values")

    class_energy = np.zeros(spatial + (classes,), dtype=np.float32)
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    for class_index, kind in enumerate(types):
        if kind == "MASK":
            energy = toads_mask_intensity_energy(
                images, modality, var_weights[:, class_index : class_index + 1], image_var
            )
        elif kind == "OBJECT":
            energy = toads_object_intensity_energy(
                images, centers[:, class_index : class_index + 1], modality,
                var_weights[:, class_index : class_index + 1], image_var,
                intensity_field=field, psf_bias=bias,
            )
        else:
            energy = toads_optimized_tissue_data_term(
                images, centers[:, class_index : class_index + 1],
                exponents[:, class_index : class_index + 1], modality,
                var_weights[:, class_index : class_index + 1], image_var,
                intensity_field=field, psf_bias=bias,
            )
        class_energy[interior + (class_index,)] = energy[interior + (0,)]

    coefficient = np.float32(prior_coefficient)
    smooth = np.float32(smoothing)
    prior_scale = np.float32(prior_factor)
    if not np.isfinite(coefficient) or not np.isfinite(smooth) or smooth < 0 or not np.isfinite(prior_scale):
        raise ValueError("membership coefficients must be finite and smoothing nonnegative")
    if coefficient > 0:
        inverse_prior, objective_prior = toads_inverse_shape_prior_energy(
            shapes, shape_prior, similarity, coefficient
        )
        relation_shapes = shapes
    else:
        inverse_prior = np.zeros_like(class_energy)
        objective_prior = np.zeros_like(class_energy)
        relation_shapes = np.zeros_like(shapes)
    factors, suppressed = toads_shape_relation_factors(relation_shapes, labels, prior_scale)
    class_energy[suppressed] = np.float32(0.0)

    updated = initial.copy()
    denominator_map = np.zeros(spatial, dtype=np.float32)
    objective = np.float32(0.0)
    threshold, reciprocal_cap = np.float32(1.0e-10), np.float32(1.0e10)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at = np.zeros(classes, dtype=np.float32)
                for class_index in range(classes):
                    if suppressed[x, y, z, class_index]:
                        updated[x, y, z, class_index] = np.float32(0.0)
                        continue
                    num = class_energy[x, y, z, class_index]
                    objective_value = num
                    if smooth > 0:
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other == class_index:
                                    continue
                                neighbor = updated[x + dx, y + dy, z + dz, other]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(neighbor * neighbor)
                                )
                        num = np.float32(
                            num + np.float32(np.float32(smooth_prior[class_index] * smooth) * neighbor_sum)
                        )
                        objective_increment = np.float32(
                            np.float32(np.float32(np.float32(0.5) * smooth_prior[class_index]) * smooth)
                            * neighbor_sum
                        )
                        objective_value = np.float32(objective_value + objective_increment)
                    num = np.float32(num + inverse_prior[x, y, z, class_index])
                    objective_value = np.float32(
                        objective_value + objective_prior[x, y, z, class_index]
                    )
                    reciprocal = np.float32(np.float32(1.0) / num) if num > threshold else reciprocal_cap
                    value = np.float32(factors[x, y, z, class_index] * reciprocal)
                    updated[x, y, z, class_index] = value
                    denominator = np.float32(denominator + value)
                    objective_at[class_index] = objective_value
                denominator_map[x, y, z] = denominator
                if denominator > 0:
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
                        contribution = np.float32(
                            objective_at[class_index] * updated[x, y, z, class_index]
                        )
                        contribution = np.float32(
                            contribution * updated[x, y, z, class_index]
                        )
                        objective = np.float32(objective + contribution)
    return ToadsAtlasTypeMembershipPass(types, class_energy, updated, denominator_map, objective)
