"""Compose raw descriptor data terms with one LesionTOADS membership pass.

This closes a fixed-input seam between the source-active descriptor energy
dispatch and source-ordered in-place membership update. Statistics, relation
map updates, spatial distance modulation, adaptive field fitting, and stopping
remain inputs or later workflow stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_active_class_energy import (
    LesionToadsActiveClassEnergies,
    compute_lesion_toads_active_class_energies,
)
from .lesion_descriptor_membership import run_lesion_toads_descriptor_membership_pass
from .lesion_objective_membership_pass import LesionObjectiveMembershipPass


@dataclass(frozen=True)
class LesionToadsDescriptorIterationSlice:
    """Outputs from raw-energy calculation followed by one membership sweep."""

    data_energies: LesionToadsActiveClassEnergies
    membership: LesionObjectiveMembershipPass
    stage_order: tuple[str, ...]


def run_lesion_toads_descriptor_iteration_slice(
    images: ArrayLike,
    initial_class_memberships: ArrayLike,
    descriptor_types: Sequence[str],
    class_centroids: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    channel_weights: ArrayLike,
    shape_values: ArrayLike,
    smoothing_priors: ArrayLike,
    shape_priors: ArrayLike,
    class_similarities: ArrayLike,
    lesion_similarities: ArrayLike,
    black_hole_similarities: ArrayLike,
    relation_types: Sequence[str],
    *,
    outlier_distance: float,
    smoothing: float,
    prior_coefficient: float,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    use_field: bool = False,
    intensity_fields: ArrayLike | None = None,
    use_lesion_weight: bool = False,
    flair_channel: int = -1,
    lesion_power_exponents: ArrayLike | None = None,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionToadsDescriptorIterationSlice:
    """Calculate active data terms then apply one fixed-input descriptor sweep.

    The function intentionally accepts already-estimated class statistics,
    shapes, similarities, and direct/prior relation decisions. Those states
    are held fixed in this bounded composition and are not estimated here.
    """
    values = np.asarray(images, dtype=np.float32)
    memberships = np.asarray(initial_class_memberships, dtype=np.float32)
    if values.ndim != 4 or memberships.ndim != 4 or values.shape[:3] != memberships.shape[:3]:
        raise ValueError("images and memberships must be matching spatial 4-D arrays")
    if memberships.shape[-1] != len(descriptor_types):
        raise ValueError("membership class count must match descriptor_types")
    energies = compute_lesion_toads_active_class_energies(
        values,
        class_centroids,
        lesion_centroid,
        black_hole_centroid,
        channel_weights,
        descriptor_types,
        outlier_distance=outlier_distance,
        use_field=use_field,
        intensity_fields=intensity_fields,
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
        lesion_power_exponents=lesion_power_exponents,
    )
    zero_map = np.zeros(values.shape[:3], dtype=np.float32)
    updated = run_lesion_toads_descriptor_membership_pass(
        memberships,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        energies.class_energy,
        energies.lesion_energy,
        energies.black_hole_energy,
        zero_map,
        shape_values,
        smoothing_priors,
        descriptor_types=descriptor_types,
        relation_types=relation_types,
        smoothing=smoothing,
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        black_hole_similarities=black_hole_similarities,
        prior_coefficient=prior_coefficient,
        zero=zero,
        infinity=infinity,
    )
    return LesionToadsDescriptorIterationSlice(
        data_energies=energies,
        membership=updated,
        stage_order=("computeActiveClassDataEnergies", "updateDescriptorMemberships", "normalizeMemberships"),
    )
