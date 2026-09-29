from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    run_lesion_toads_lesion_smoothing_sweep,
    update_lesion_toads_smoothing_energies,
)


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_lesion_sweep/java_reference.json"
PRIOR_REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_lesion_prior_sweep/java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_TYPES = ("OBJECT", "LESION", "OBJECT")
RELATIONS = ("direct", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)
FACES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    class_terms = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        class_terms[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    lesion_terms = (0.75 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    shapes[2, 2, 2, 1] = 0.0
    return memberships, lesion, class_terms, lesion_terms, shapes


def _reverse_order_control(
    initial: np.ndarray,
    lesion_input: np.ndarray,
    class_data: np.ndarray,
    lesion_data: np.ndarray,
    shapes: np.ndarray,
    smoothing: float,
) -> tuple[np.ndarray, np.ndarray]:
    current, lesion = initial.copy(), lesion_input.copy()
    coords = list(itertools.product(*(range(1, size - 1) for size in current.shape[:3])))
    zero_faces = np.zeros(6, dtype=np.float32)
    for x, y, z in reversed(coords):
        points = [(x + dx, y + dy, z + dz) for dx, dy, dz in FACES]
        class_faces = np.stack([current[point] for point in points])
        lesion_faces = np.asarray([lesion[point] for point in points], dtype=np.float32)
        denominator = np.float32(0)
        for k, kind in enumerate(CLASS_TYPES):
            if shapes[x, y, z, k] <= 0:
                current[x, y, z, k] = 0
                if kind == "LESION":
                    lesion[x, y, z] = 0
            else:
                update = update_lesion_toads_smoothing_energies(
                    class_faces, lesion_faces, zero_faces,
                    class_index=k, class_type=kind,
                    smoothing_prior=float(SMOOTHING_PRIORS[k]), smoothing=smoothing,
                    class_energy=0, lesion_energy=0, black_hole_energy=0,
                    class_membership_data_term=float(class_data[x, y, z, k]),
                    lesion_membership_data_term=(float(lesion_data[x, y, z]) if kind == "LESION" else 0),
                    black_hole_membership_data_term=0, outmask_membership_data_term=0,
                )
                current[x, y, z, k] = convert_lesion_toads_class_membership(
                    float(update.class_membership_data_term), float(shapes[x, y, z, k]),
                    relation_type=RELATIONS[k],
                )
                if kind == "LESION":
                    lesion[x, y, z] = convert_lesion_toads_class_membership(
                        float(update.lesion_membership_data_term), float(shapes[x, y, z, k]),
                        relation_type=RELATIONS[k],
                    )
            if kind == "LESION":
                denominator = np.float32(denominator + lesion[x, y, z])
            denominator = np.float32(denominator + current[x, y, z, k])
        if denominator > 0:
            for k, kind in enumerate(CLASS_TYPES):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if kind == "LESION":
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
    return current, lesion


def test_object_lesion_sweep_matches_selected_java_in_place_references() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    memberships, lesion, class_data, lesion_data, shapes = _inputs(
        tuple(reference["shape"]), reference["classes"]
    )
    initial_memberships, initial_lesion = memberships.copy(), lesion.copy()
    result = run_lesion_toads_lesion_smoothing_sweep(
        memberships,
        lesion,
        class_data,
        lesion_data,
        shapes,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
    border = np.ones(memberships.shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], initial_memberships[border])
    np.testing.assert_array_equal(result.lesion_membership[border], initial_lesion[border])
    assert result.lesion_membership[2, 2, 2] == 0
    reverse = _reverse_order_control(
        initial_memberships,
        initial_lesion,
        class_data,
        lesion_data,
        shapes,
        float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )
    assert np.any(result.class_memberships[~border] != reverse[0][~border])
    np.testing.assert_array_equal(memberships, initial_memberships)
    np.testing.assert_array_equal(lesion, initial_lesion)


def test_lesion_sweep_rejects_multiple_lesion_classes() -> None:
    memberships, lesion, class_data, lesion_data, shapes = _inputs((3, 3, 3), 3)
    with pytest.raises(ValueError, match="exactly one LESION class"):
        run_lesion_toads_lesion_smoothing_sweep(
            memberships,
            lesion,
            class_data,
            lesion_data,
            shapes,
            SMOOTHING_PRIORS,
            ("OBJECT", "LESION", "LESION"),
            RELATIONS,
            smoothing=0.2,
        )


def test_lesion_prior_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(PRIOR_REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    memberships, lesion, class_data, lesion_data, shapes = _inputs(
        tuple(reference["shape"]), reference["classes"]
    )
    shape_priors = np.asarray([0.25, 0.375, 0.5], dtype=np.float32)
    class_similarities = np.asarray(
        [[0.0, 0.125, 0.25], [0.5, 0.0, 0.375], [0.25, 0.625, 0.0]],
        dtype=np.float32,
    )
    lesion_similarities = np.asarray([0.25, 0.0, 0.5], dtype=np.float32)
    result = run_lesion_toads_lesion_smoothing_sweep(
        memberships,
        lesion,
        class_data,
        lesion_data,
        shapes,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        prior_coefficient=float(
            np.asarray(reference["prior_coefficient_bits"], dtype=np.uint32).view(np.float32)
        ),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]

    no_prior = run_lesion_toads_lesion_smoothing_sweep(
        memberships,
        lesion,
        class_data,
        lesion_data,
        shapes,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=0.2,
    )
    explicit_zero_prior = run_lesion_toads_lesion_smoothing_sweep(
        memberships,
        lesion,
        class_data,
        lesion_data,
        shapes,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=0.2,
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        prior_coefficient=0.0,
    )
    np.testing.assert_array_equal(no_prior.class_memberships, explicit_zero_prior.class_memberships)
    np.testing.assert_array_equal(no_prior.lesion_membership, explicit_zero_prior.lesion_membership)

    with pytest.raises(ValueError, match="required for a positive prior_coefficient"):
        run_lesion_toads_lesion_smoothing_sweep(
            memberships,
            lesion,
            class_data,
            lesion_data,
            shapes,
            SMOOTHING_PRIORS,
            CLASS_TYPES,
            RELATIONS,
            smoothing=0.2,
            prior_coefficient=0.3,
        )
