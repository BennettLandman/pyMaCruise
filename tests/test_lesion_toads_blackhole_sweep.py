from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_blackhole_smoothing_sweep


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_blackhole_sweep/java_reference.json"
PRIOR_REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_blackhole_prior_sweep/java_reference.json"
DISTANCE_REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_blackhole_distance_sweep/java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_TYPES = ("OBJECT", "BLACKHOLE", "OBJECT")
RELATIONS = ("direct", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    class_data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        class_data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    black_hole = (0.5 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    lesion_data = (0.75 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    black_hole_data = (0.5 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    shapes[2, 2, 2, 1] = 0.0
    return memberships, lesion, black_hole, class_data, lesion_data, black_hole_data, shapes


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_blackhole_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    before = tuple(value.copy() for value in inputs[:3])
    result = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(result.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]
    border = np.ones(inputs[0].shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], before[0][border])
    np.testing.assert_array_equal(result.lesion_membership[border], before[1][border])
    np.testing.assert_array_equal(result.black_hole_membership[border], before[2][border])
    np.testing.assert_array_equal(inputs[0], before[0])
    np.testing.assert_array_equal(inputs[1], before[1])
    np.testing.assert_array_equal(inputs[2], before[2])


def test_blackhole_sweep_rejects_multiple_blackhole_classes() -> None:
    inputs = _inputs((3, 3, 3), 3)
    with pytest.raises(ValueError, match="exactly one BLACKHOLE class"):
        run_lesion_toads_blackhole_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            ("OBJECT", "BLACKHOLE", "BLACKHOLE"),
            RELATIONS,
            smoothing=0.2,
        )


def test_blackhole_prior_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(PRIOR_REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    shape_priors = np.asarray([0.25, 0.375, 0.5], dtype=np.float32)
    class_similarities = np.asarray(
        [[0.0, 0.125, 0.25], [0.5, 0.0, 0.375], [0.25, 0.625, 0.0]],
        dtype=np.float32,
    )
    lesion_similarities = np.asarray([0.25, 0.0, 0.5], dtype=np.float32)
    black_hole_similarities = np.asarray([0.125, 0.0, 0.625], dtype=np.float32)
    result = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        black_hole_similarities=black_hole_similarities,
        prior_coefficient=float(
            np.asarray(reference["prior_coefficient_bits"], dtype=np.uint32).view(np.float32)
        ),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(result.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]

    no_prior = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs, SMOOTHING_PRIORS, CLASS_TYPES, RELATIONS, smoothing=0.2
    )
    explicit_zero_prior = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=0.2,
        shape_priors=shape_priors,
        class_similarities=class_similarities,
        lesion_similarities=lesion_similarities,
        black_hole_similarities=black_hole_similarities,
        prior_coefficient=0.0,
    )
    np.testing.assert_array_equal(no_prior.class_memberships, explicit_zero_prior.class_memberships)
    np.testing.assert_array_equal(no_prior.lesion_membership, explicit_zero_prior.lesion_membership)
    np.testing.assert_array_equal(no_prior.black_hole_membership, explicit_zero_prior.black_hole_membership)

    with pytest.raises(ValueError, match="required for a positive prior_coefficient"):
        run_lesion_toads_blackhole_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            CLASS_TYPES,
            RELATIONS,
            smoothing=0.2,
            prior_coefficient=0.3,
        )


def test_blackhole_distance_sweep_matches_selected_java_asymmetry() -> None:
    reference = json.loads(DISTANCE_REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    shape = tuple(reference["shape"])
    inputs = _inputs(shape, reference["classes"])
    memberships, lesion, black_hole, class_data, lesion_data, black_hole_data, shapes = inputs
    vent = np.full(shape, -1, dtype=np.int8)
    gm = np.full(shape, -1, dtype=np.int8)
    brainstem = np.full(shape, -1, dtype=np.int8)
    for x, y, z in np.ndindex(shape):
        coordinate_sum = x + y + z
        if coordinate_sum % 3 == 0:
            vent[x, y, z] = (x + y) % 4
        if coordinate_sum % 3 != 0 and (x + 2 * y + z) % 2 == 0:
            gm[x, y, z] = coordinate_sum % 5
        if (x + 2 * y + z) % 2 == 0:
            brainstem[x, y, z] = coordinate_sum % 6
    before = tuple(value.copy() for value in inputs[:3])
    results = {}
    for expected_case in reference["cases"]:
        gm_flag = expected_case["name"] == "gm_enabled"
        result = run_lesion_toads_blackhole_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            CLASS_TYPES,
            RELATIONS,
            smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
            ventricle_distances=vent,
            gray_matter_distances=gm,
            brainstem_distances=brainstem,
            max_ventricle_distance=3,
            max_gray_matter_distance=4,
            max_brainstem_distance=5,
            gm_flag=gm_flag,
            spread=1.25,
        )
        results[expected_case["name"]] = result
        for expected in expected_case["voxels"]:
            point = tuple(expected["coordinate"])
            assert _bits(result.class_memberships[point]) == expected["class_bits"]
            assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
            assert int(np.asarray(result.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]

    no_distance = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs, SMOOTHING_PRIORS, CLASS_TYPES, RELATIONS, smoothing=0.2
    )
    baseline = json.loads(REFERENCE.read_text())
    for expected in baseline["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(no_distance.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(no_distance.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(no_distance.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    assert np.any(
        results["gm_enabled"].black_hole_membership[interior]
        != results["gm_disabled"].black_hole_membership[interior]
    )
    for original, after in zip(before, inputs[:3], strict=True):
        np.testing.assert_array_equal(original, after)
