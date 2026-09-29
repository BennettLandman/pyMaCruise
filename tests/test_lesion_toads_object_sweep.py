from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    run_lesion_toads_object_membership_sweep,
    update_lesion_toads_smoothing_energies,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_object_sweep/java_reference.json"
)
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
FACES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _synthetic_inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    energies = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        energies[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
    return memberships, energies


def _run_reverse_order(
    initial: np.ndarray,
    energies: np.ndarray,
    priors: np.ndarray,
    smoothing: float,
) -> np.ndarray:
    """Sensitivity control: use same scalar kernel, reverse only scan order."""
    current = initial.copy()
    shape = current.shape[:3]
    coords = list(itertools.product(*(range(1, size - 1) for size in shape)))
    for x, y, z in reversed(coords):
        points = [(x + dx, y + dy, z + dz) for dx, dy, dz in FACES]
        class_faces = np.stack([current[point] for point in points])
        zero_faces = np.zeros(6, dtype=np.float32)
        denominator = np.float32(0.0)
        for k in range(current.shape[3]):
            update = update_lesion_toads_smoothing_energies(
                class_faces,
                zero_faces,
                zero_faces,
                class_index=k,
                class_type="OBJECT",
                smoothing_prior=float(priors[k]),
                smoothing=smoothing,
                class_energy=0,
                lesion_energy=0,
                black_hole_energy=0,
                class_membership_data_term=float(energies[x, y, z, k]),
                lesion_membership_data_term=0,
                black_hole_membership_data_term=0,
                outmask_membership_data_term=0,
            )
            current[x, y, z, k] = convert_lesion_toads_class_membership(
                float(update.class_membership_data_term), 1.0, relation_type="direct"
            )
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > 0:
            current[x, y, z] = np.asarray(current[x, y, z] / denominator, dtype=np.float32)
    return current


def test_direct_object_sweep_matches_java_in_place_voxel_references() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    priors = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)
    initial, energies = _synthetic_inputs(tuple(reference["shape"]), reference["classes"])
    unchanged = initial.copy()
    actual = run_lesion_toads_object_membership_sweep(
        initial, energies, priors,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )

    for coordinate, expected in zip(reference["voxels"], reference["voxels"], strict=True):
        point = tuple(coordinate["coordinate"])
        assert _bits(actual[point]) == expected["class_bits"]
    boundary = np.ones(initial.shape[:3], dtype=bool)
    boundary[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(actual[boundary], unchanged[boundary])
    assert np.any(actual[~boundary] != _run_reverse_order(unchanged, energies, priors, 0.2)[~boundary])
    np.testing.assert_array_equal(initial, unchanged)


@pytest.mark.parametrize(
    ("initial", "energies", "priors", "smoothing", "message"),
    [
        (np.zeros((2, 3, 3, 1)), np.zeros((2, 3, 3, 1)), np.ones(1), 0.2, "spatial sizes >= 3"),
        (np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 2)), np.ones(1), 0.2, "must match"),
        (np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.ones(1), 0.0, "must be finite and positive"),
    ],
)
def test_object_sweep_rejects_unsupported_or_misaligned_inputs(
    initial: np.ndarray,
    energies: np.ndarray,
    priors: np.ndarray,
    smoothing: float,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        run_lesion_toads_object_membership_sweep(initial, energies, priors, smoothing=smoothing)
