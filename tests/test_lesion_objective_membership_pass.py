from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_objective_membership_pass import run_lesion_toads_objective_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_objective_membership_pass"
REFERENCE = ROOT / "java_reference.json"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
TYPES = ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK")
RELATIONS = ("direct", "prior", "direct", "prior")


def _expected(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    unsigned = np.asarray([int(value) & 0xFFFFFFFF for value in bits], dtype=np.uint32)
    return unsigned.view(np.float32).reshape(shape)


def test_supported_membership_maps_match_fixed_input_java_harness() -> None:
    image = nib.load(str(ROOT / "inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    inputs = np.asarray(image.dataobj, dtype=np.float32)
    assert inputs.shape == (5, 4, 5, 18)
    input_sha = hashlib.sha256(inputs.astype("<f4", copy=False).view("<u4").tobytes()).hexdigest()

    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["input_sha256"] == input_sha
    assert reference["shape"] == [5, 4, 5]
    assert reference["class_types"] == list(TYPES)

    result = run_lesion_toads_objective_membership_pass(
        inputs[..., 0:4],
        inputs[..., 4],
        inputs[..., 5],
        inputs[..., 6],
        inputs[..., 7:11],
        inputs[..., 11],
        inputs[..., 12],
        inputs[..., 13],
        inputs[..., 14:18],
        np.asarray([0.375, 0.5, 0.625, 0.75], dtype=np.float32),
        TYPES,
        RELATIONS,
        smoothing=0.25,
        shape_priors=np.asarray([0.125, 0.25, 0.375, 0.5], dtype=np.float32),
        class_similarities=np.asarray(
            [
                [0.0, 0.3125, 0.4375, 0.5625],
                [0.25, 0.0, 0.5, 0.625],
                [0.375, 0.3125, 0.0, 0.6875],
                [0.1875, 0.4375, 0.5625, 0.0],
            ],
            dtype=np.float32,
        ),
        lesion_similarities=np.asarray([0.6875, 0.5625, 0.4375, 0.3125], dtype=np.float32),
        black_hole_similarities=np.asarray([0.25, 0.375, 0.5, 0.625], dtype=np.float32),
        prior_coefficient=0.375,
    )

    actual_arrays = {
        "class_memberships_bits": result.class_memberships,
        "lesion_membership_bits": result.lesion_membership,
        "black_hole_membership_bits": result.black_hole_membership,
        "outlier_membership_bits": result.outlier_membership,
        "denominator_bits": result.denominator,
    }
    for name, actual in actual_arrays.items():
        expected = _expected(reference[name], actual.shape)
        np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32), err_msg=name)

    # The source scans only the six-face interior; all boundaries retain inputs
    # and generated maps stay zero there.
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    for axis in range(3):
        for edge in (0, -1):
            key = [slice(None)] * 3
            key[axis] = edge
            key = tuple(key)
            np.testing.assert_array_equal(result.class_memberships[key], inputs[..., 0:4][key])
            np.testing.assert_array_equal(result.lesion_membership[key], inputs[..., 4][key])
            np.testing.assert_array_equal(result.black_hole_membership[key], inputs[..., 5][key])
            np.testing.assert_array_equal(result.outlier_membership[key], inputs[..., 6][key])
            assert not np.any(result.denominator[key])
    assert np.all(result.denominator[interior] > 0)


def test_objective_membership_pass_rejects_shape_skipped_interior_classes() -> None:
    inputs = np.asarray(nib.load(str(ROOT / "inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    shapes = inputs[..., 14:18].copy()
    shapes[2, 1, 2, 0] = 0.0
    with np.testing.assert_raises_regex(ValueError, "all interior classes must be active"):
        run_lesion_toads_objective_membership_pass(
            inputs[..., 0:4], inputs[..., 4], inputs[..., 5], inputs[..., 6],
            inputs[..., 7:11], inputs[..., 11], inputs[..., 12], inputs[..., 13],
            shapes, np.asarray([0.375, 0.5, 0.625, 0.75], dtype=np.float32), TYPES, RELATIONS,
            smoothing=0.25,
            shape_priors=np.asarray([0.125, 0.25, 0.375, 0.5], dtype=np.float32),
            class_similarities=np.eye(4, dtype=np.float32),
            lesion_similarities=np.ones(4, dtype=np.float32),
            black_hole_similarities=np.ones(4, dtype=np.float32),
            prior_coefficient=0.0,
        )


def test_source_faithful_mode_retains_last_energy_after_zero_shape_skip() -> None:
    inputs = np.asarray(nib.load(str(ROOT / "inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    shapes = inputs[..., 14:18].copy()
    coordinate = (2, 1, 2)
    class_index = 0

    def run(shape_values: np.ndarray):
        return run_lesion_toads_objective_membership_pass(
            inputs[..., 0:4], inputs[..., 4], inputs[..., 5], inputs[..., 6],
            inputs[..., 7:11], inputs[..., 11], inputs[..., 12], inputs[..., 13],
            shape_values, np.asarray([0.375, 0.5, 0.625, 0.75], dtype=np.float32), TYPES, RELATIONS,
            smoothing=0.25,
            shape_priors=np.asarray([0.125, 0.25, 0.375, 0.5], dtype=np.float32),
            class_similarities=np.eye(4, dtype=np.float32),
            lesion_similarities=np.ones(4, dtype=np.float32),
            black_hole_similarities=np.ones(4, dtype=np.float32),
            prior_coefficient=0.0,
            source_faithful=True,
        )

    active_result = run(shapes)
    skipped_shapes = shapes.copy()
    skipped_shapes[coordinate + (class_index,)] = 0.0
    skipped_result = run(skipped_shapes)

    assert skipped_result.class_memberships[coordinate + (class_index,)] == 0.0
    previous_coordinate = (coordinate[0], coordinate[1], coordinate[2] - 1)
    assert skipped_result.class_energy[coordinate + (class_index,)] == (
        skipped_result.class_energy[previous_coordinate + (class_index,)]
    )
    # Skipping the class leaves the method-local energy value untouched,
    # unlike evaluating the positive-shape class at that same voxel.
    assert skipped_result.class_energy[coordinate + (class_index,)] != (
        active_result.class_energy[coordinate + (class_index,)]
    )
