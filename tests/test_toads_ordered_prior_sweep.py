from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_math import (
    toads_inverse_shape_prior_energy,
    toads_ordered_membership_sweep,
    toads_shape_relation_factors,
)


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_ordered_prior_sweep" / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
SMOOTHING = np.float32(0.017)
PRIOR_FACTOR = np.float32(0.75)
SMOOTHING_PRIORS = np.asarray([0.5, 0.8, 1.1, 0.25], dtype=np.float32)


def _load(directory: str, name: str) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    energy = _load("toads_direct_membership", "data_energy")
    initial = _load("toads_direct_membership", "initial_memberships")
    shapes = _load("toads_shape_relations", "shape_priors")
    labels = _load("toads_shape_relations", "closest_labels").astype(np.int16)
    return energy, initial, shapes, labels


def _prior_terms(reference: dict, shapes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    prior_volume = toads_inverse_shape_prior_energy(
        shapes,
        reference["atlas_shape_priors"],
        reference["similarity"],
        reference["prior_coefficient"],
    )
    return prior_volume


def test_ordered_sweep_with_inverse_shape_prior_matches_java_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["shape"] == [7, 6, 5, 4]
    energy, initial, shapes, labels = _inputs()
    prior_inverse, prior_objective = _prior_terms(reference, shapes)
    updated, objective = toads_ordered_membership_sweep(
        energy,
        initial,
        shapes,
        labels,
        SMOOTHING_PRIORS,
        SMOOTHING,
        PRIOR_FACTOR,
        prior_inverse_energy=prior_inverse,
        prior_objective_energy=prior_objective,
    )
    np.testing.assert_array_equal(
        updated, np.asarray(reference["updated"], dtype=np.float32).reshape(updated.shape)
    )
    assert objective == np.float32(reference["objective"])
    np.testing.assert_array_equal(updated[0], initial[0])
    np.testing.assert_array_equal(updated[:, -1], initial[:, -1])


def test_suppressed_prior_only_class_skips_prior_objective() -> None:
    reference = json.loads(REFERENCE.read_text())
    energy, initial, shapes, labels = _inputs()
    factors, suppressed = toads_shape_relation_factors(shapes, labels, PRIOR_FACTOR)
    locations = np.argwhere(suppressed)
    assert locations.size
    x, y, z, class_index = (int(v) for v in locations[0])
    inverse, prior_objective = _prior_terms(reference, shapes)
    baseline, baseline_objective = toads_ordered_membership_sweep(
        energy, initial, shapes, labels, SMOOTHING_PRIORS, SMOOTHING, PRIOR_FACTOR,
        prior_inverse_energy=inverse, prior_objective_energy=prior_objective,
    )
    inverse[x, y, z, class_index] = np.float32(10.0)
    prior_objective[x, y, z, class_index] = np.float32(5.0)
    updated, objective = toads_ordered_membership_sweep(
        energy, initial, shapes, labels, SMOOTHING_PRIORS, SMOOTHING, PRIOR_FACTOR,
        prior_inverse_energy=inverse, prior_objective_energy=prior_objective,
    )
    np.testing.assert_array_equal(updated, baseline)
    assert objective == baseline_objective


def test_zero_prior_energy_recovers_pre_prior_sweep() -> None:
    energy, initial, shapes, labels = _inputs()
    zeros = np.zeros_like(energy)
    without, objective_without = toads_ordered_membership_sweep(
        energy, initial, shapes, labels, SMOOTHING_PRIORS, SMOOTHING, PRIOR_FACTOR
    )
    with_zeros, objective_with = toads_ordered_membership_sweep(
        energy,
        initial,
        shapes,
        labels,
        SMOOTHING_PRIORS,
        SMOOTHING,
        PRIOR_FACTOR,
        prior_inverse_energy=zeros,
        prior_objective_energy=zeros,
    )
    np.testing.assert_array_equal(with_zeros, without)
    assert objective_with == objective_without
