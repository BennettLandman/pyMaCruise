"""Selected-Java comparison for the default 11-class approximate pyramid."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_approximate_initial_membership import (
    initialize_toads_default_approximate_initial_memberships,
)
from pymacruise.toads_default_approximate_registration import (
    register_toads_default_approximate_atlas_pyramid,
)
from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_approximate_registration"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultApproximateRegistrationReferenceHarness.java"
HARNESS_SHA256 = "973592718d4d55e20aed157909e5deb566c5cf4a3a552e76b6eddaf1a11abd0e"
DESCRIPTOR = Path(__file__).parents[1] / "data" / "atlases" / "cruise-atlas-10obj-toads2008.txt"


def _bits(value: np.ndarray | np.float32) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def _decode_java_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32)


def _ordered_float_bits(values: np.ndarray | np.float32) -> np.ndarray:
    signed = np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(np.int64)
    return np.where(signed < 0, -signed, signed + (1 << 31))


def _max_ulp_difference(actual: np.ndarray | np.float32, java_bits: list[int]) -> int:
    expected = _decode_java_bits(java_bits)
    return int(np.max(np.abs(_ordered_float_bits(actual) - _ordered_float_bits(expected))))


def _load_inputs():
    images_nii = nib.load(str(ROOT / "images.nii.gz"))
    atlas_nii = nib.load(str(ROOT / "atlas_shapes.nii.gz"))
    for image in (images_nii, atlas_nii):
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(
            nib.affines.voxel_sizes(image.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6,
        )
    assert images_nii.shape == (9, 8, 7, 2)
    assert atlas_nii.shape == (11, 10, 9, 11)
    assert not np.array_equal(images_nii.shape[:3], atlas_nii.shape[:3])
    return tuple(np.asarray(image.dataobj, dtype=np.float32) for image in (images_nii, atlas_nii))


def _centroids() -> np.ndarray:
    result = np.asarray(
        [[np.float32(.19) + np.float32(.021) * k + np.float32(.013) * c for k in range(11)] for c in range(2)],
        dtype=np.float32,
    )
    result[0, 2] = result[0, 1]
    return result


def test_default_approximate_pyramid_matches_selected_java_reference() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert JAVA["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert JAVA["preprocess_source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert JAVA["class_types"] == list(DEFAULT_ATLAS_CLASS_TYPES)
    assert JAVA["call_order"] == [
        "initialMemberships", "setImageInfo", "initShapeRegistration",
        "setTransform(preprocess)", "updateShapeRegistration", "registerShapesPyramid",
    ]
    assert JAVA["pyramid_dispatch"] == "registerShapesApproxPyramid"
    assert JAVA["scale_order"] == [4, 2]
    assert JAVA["registered_class_indices"] == list(range(1, 11))
    assert hashlib.sha256(DESCRIPTOR.read_bytes()).hexdigest() == JAVA["descriptor_sha256"]

    images, atlas_shapes = _load_inputs()
    originals = tuple(values.copy() for values in (images, atlas_shapes))
    init = initialize_toads_default_approximate_initial_memberships(
        images, _centroids(), approx=JAVA["approx"],
    )
    assert _bits(init.memberships) == JAVA["initial_memberships_bits"]
    np.testing.assert_array_equal(
        init.best_class_indices.ravel().astype(int), JAVA["initial_best_indices"],
    )
    result = register_toads_default_approximate_atlas_pyramid(
        init,
        atlas_shapes,
        image_resolutions_mm=(.9, 1.2, 1.45),
        template_resolutions_mm=(.9, 1.2, 1.45),
        initial_transform=(.004, -.003, .002, .04, -.03, .02),
        levels=JAVA["levels"],
        max_iterations=3,
    )
    assert len(result.levels) == 2
    for actual, reference in zip(result.levels, JAVA["stages"]):
        assert actual.scale == reference["scale"]
        assert list(actual.memberships.shape[:3]) == reference["dims"]
        assert _bits(actual.memberships) == reference["memberships_bits"]
        np.testing.assert_array_equal(
            actual.best_class_indices.ravel().astype(int), reference["best_indices"],
        )
        assert _bits(actual.initial_transform) == reference["initial_transform_bits"]
        assert len(actual.iteration_coefficients) == len(reference["iterations"])
        assert len(actual.iteration_trials) == len(reference["iterations"])
        for iteration_index, (coefficients, trials, java_iteration, py_iteration) in enumerate(zip(
            actual.iteration_coefficients, actual.iteration_trials,
            reference["iterations"], actual.result.iterations, strict=True,
        )):
            if actual.scale == 4:
                assert _bits(coefficients[0]) == java_iteration["gradient"]
                assert _bits(coefficients[1]) == java_iteration["hessian"]
            elif actual.scale == 2 and iteration_index == 0:
                # Only the first scale-2 coefficient pass shows bounded
                # float32 differences against Java's scalar interpolation.
                assert _max_ulp_difference(coefficients[0], java_iteration["gradient"]) == 23
                assert _max_ulp_difference(coefficients[1], java_iteration["hessian"]) == 1
            else:
                assert _bits(coefficients[0]) == java_iteration["gradient"]
                assert _bits(coefficients[1]) == java_iteration["hessian"]
            assert len(trials) == len(java_iteration["trials"])
            for trial_index, ((trial_transform, trial_energy), java_trial) in enumerate(
                zip(trials, java_iteration["trials"])
            ):
                if actual.scale == 4:
                    assert _bits(trial_transform) == java_trial["transform"]
                    assert _bits(trial_energy) == [java_trial["energy"]]
                elif actual.scale == 2:
                    trial_delta = np.abs(
                        trial_transform - _decode_java_bits(java_trial["transform"])
                    )
                    # The rejected lambda=1 scale-2 trial has the largest
                    # parameter drift; later trials inherit the accepted-state
                    # drift but remain within this bound.
                    assert float(np.max(trial_delta)) <= 7.629395e-6
                    if iteration_index == 0 and trial_index == 0:
                        assert _max_ulp_difference(trial_energy, [java_trial["energy"]]) == 3
                    else:
                        assert _bits(trial_energy) == [java_trial["energy"]]
                else:
                    assert _bits(trial_transform) == java_trial["transform"]
                    assert _bits(trial_energy) == [java_trial["energy"]]
            assert _bits(py_iteration.lambda_value) == [java_iteration["lambda"]]
            assert _bits(py_iteration.energy) == [java_iteration["energy"]]
            assert _bits(py_iteration.relative_change) == [java_iteration["relative_change"]]
            assert py_iteration.accepted == java_iteration["accepted"]
            assert py_iteration.rejected_trials == java_iteration["rejected_trials"]
            if actual.scale == 4:
                assert _bits(py_iteration.transform) == java_iteration["transform"]
            elif actual.scale == 2:
                # The accumulated parameter vector retains the small first
                # scale-2 gradient difference; scalar state/control remains exact.
                assert float(np.max(np.abs(
                    py_iteration.transform - _decode_java_bits(java_iteration["transform"])
                ))) <= 3.814698e-6
            assert _bits(py_iteration.transform) == _bits(trials[-1][0])
        assert _bits(actual.result.iterations[-1].lambda_value) == [reference["lambda_bits"]]
        assert actual.result.stop_reason == reference["stop_reason"]
        assert _bits(actual.result.transform) == _bits(actual.result.iterations[-1].transform)
        if actual.scale == 4:
            assert _bits(actual.result.transform) == reference["transform_bits"]
            assert _bits(actual.result.energy) == [reference["energy_bits"]]
        else:
            assert float(np.max(np.abs(
                actual.result.transform - _decode_java_bits(reference["transform_bits"])
            ))) <= 3.814698e-6
            assert _bits(actual.result.energy) == [reference["energy_bits"]]
    assert float(np.max(np.abs(
        result.transform - _decode_java_bits(JAVA["final_transform_bits"])
    ))) <= 3.814698e-6
    assert _bits(result.transform) == _bits(result.levels[-1].result.transform)
    assert len(result.shape_bounds) == 11
    assert any(bound.support_count for bound in result.shape_bounds[1:])
    assert result.rotation.shape == (3, 3)
    for actual, original in zip((images, atlas_shapes), originals):
        np.testing.assert_array_equal(actual, original)
