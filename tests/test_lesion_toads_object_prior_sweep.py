from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    run_lesion_toads_object_membership_sweep,
    run_lesion_toads_object_prior_membership_sweep,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_object_prior_sweep/java_reference.json"
)
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    shapes[2, 2, 2, 1] = 0.0
    shape_priors = np.asarray([0.25, 0.375, 0.5], dtype=np.float32)
    similarities = np.asarray(
        [[0.0, 0.125, 0.25], [0.5, 0.0, 0.375], [0.25, 0.625, 0.0]],
        dtype=np.float32,
    )
    smoothing_priors = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)
    relations = ("direct", "prior", "prior")
    return memberships, data, shapes, shape_priors, similarities, smoothing_priors, relations


def test_object_prior_sweep_matches_selected_java_in_place_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    result = run_lesion_toads_object_prior_membership_sweep(
        *inputs,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
        prior_coefficient=float(np.asarray(reference["prior_coefficient_bits"], dtype=np.uint32).view(np.float32)),
    )
    for expected in reference["voxels"]:
        coordinate = tuple(expected["coordinate"])
        assert _bits(result[coordinate]) == expected["class_bits"]


def test_zero_prior_all_direct_matches_previous_object_sweep() -> None:
    inputs = _inputs((5, 4, 6), 3)
    memberships, data, shapes, shape_priors, similarities, smoothing_priors, _ = inputs
    shapes.fill(1.0)
    direct = ("direct", "direct", "direct")
    with_prior_adapter = run_lesion_toads_object_prior_membership_sweep(
        memberships,
        data,
        shapes,
        shape_priors,
        similarities,
        smoothing_priors,
        direct,
        smoothing=0.2,
        prior_coefficient=0.0,
    )
    previous = run_lesion_toads_object_membership_sweep(
        memberships, data, smoothing_priors, smoothing=0.2
    )
    np.testing.assert_array_equal(with_prior_adapter, previous)
