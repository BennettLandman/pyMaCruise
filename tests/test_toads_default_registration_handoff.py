from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_registration_handoff import (
    toads_default_atlas_registration_handoff,
)
from pymacruise.toads_registration_sampling import register_toads_atlas_rigid
from pymacruise.toads_shapes import compute_toad_shape_bounds

from test_toads_default_atlas_dispatch import _fixture, _run


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_registration_handoff"
JAVA = json.loads((ROOT / "java_reference.json").read_text())


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _array_from_bits(values: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape)


def _arguments():
    data, parameters, cases, reference = _fixture()
    source = nib.load(str(ROOT / "atlas_shapes.nii.gz"))
    assert nib.aff2axcodes(source.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(source.affine), (0.9, 1.2, 1.45), rtol=0, atol=1e-6)
    subject = nib.load(str(Path(__file__).parent / "data" / "synthetic" / "toads_default_atlas_dispatch" / "images.nii.gz"))
    assert nib.aff2axcodes(subject.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(subject.affine), (0.9, 1.2, 1.45), rtol=0, atol=1e-6)
    np.testing.assert_allclose(
        source.affine @ np.asarray([1.0, 1.0, 1.0, 1.0]),
        subject.affine @ np.asarray([0.0, 0.0, 0.0, 1.0]),
        rtol=0,
        atol=2e-6,
    )
    atlas_shapes = np.asarray(source.dataobj, dtype=np.float32)
    assert atlas_shapes.shape == (11, 10, 9, 11)
    assert atlas_shapes.shape[:3] != data["images"].shape[:3]
    case = next(
        item for item in cases
        if item["field"] and item["psf"] and item["smoothing"] and item["prior"]
    )
    options = {
        "image_values": data["images"],
        "initial_memberships": data["initial_memberships"],
        "centroids": parameters["centroids"],
        "exponent_by_channel_class": parameters["exponent_by_channel_class"],
        "modality_weights": parameters["modality_weights"],
        "class_variance_weights": parameters["class_variance_weights"],
        "image_variances": parameters["image_variances"],
        "transformed_shapes": data["transformed_shapes"],
        "closest_labels": data["closest_labels"],
        "smoothing_priors": parameters["smoothing_priors"],
        "centroid_similarity": parameters["centroid_similarity"],
        "atlas_shape_priors": parameters["atlas_shape_priors"],
        "smoothing": 0.043,
        "prior_factor": parameters["prior_factor"],
        "prior_coefficient": 0.16,
        "atlas_shapes": atlas_shapes,
        "image_resolutions_mm": (0.9, 1.2, 1.45),
        "template_resolutions_mm": (0.9, 1.2, 1.45),
        "initial_transform": (0.004, -0.003, 0.002, 0.04, -0.03, 0.02),
        "intensity_field": data["field"],
        "psf_bias": data["psf_bias"],
        # Four steps distinguish the freshly updated membership registration
        # trajectory from the initial-membership baseline on this phantom.
        "max_iterations": 4,
        "energy_tolerance": 0.0,
        "minimum_lambda": 1e-3,
    }
    return data, parameters, reference, case, options, subject.affine


def test_default_membership_outputs_feed_rigid_registration_in_wrapper_order() -> None:
    data, parameters, membership_reference, case, options, subject_affine = _arguments()
    assert JAVA["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert JAVA["atlas_source_sha256"] == "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db"
    assert JAVA["call_order"][:5] == [
        "computeExactMemberships", "energyRatio", "registerShapes",
        "convergenceCheck", "iterationLimitCheck",
    ]
    assert JAVA["align_shapes_default"] is True
    assert JAVA["registration_dispatch"] == "SINGLE/registerAllShapes"
    assert JAVA["registration_resets"] == ["lambda=INIT_LAMBDA", "Nturn=0", "itPlus=0"]
    assert JAVA["registration_stop"] == ["n < itMax", "diff > minEdiff", "lambda > minLambda"]
    assert JAVA["registration_finalize"] == ["computeRotation(transform)", "computeTransformedShapeBoundingBox()"]

    expected_membership = _run(data, parameters, case)
    assert _bits(expected_membership.memberships) == case["membership_bits"]
    assert _bits(expected_membership.class_data_energy) == case["class_energy_bits"]
    assert _bits(expected_membership.denominator) == case["denominator_bits"]
    expected_membership_from_java = _array_from_bits(
        case["membership_bits"], data["initial_memberships"].shape
    )

    before = {key: value.copy() for key, value in data.items()}
    statistical_inputs_before = {
        key: np.asarray(value).copy()
        for key, value in options.items()
        if isinstance(value, np.ndarray)
    }
    result = toads_default_atlas_registration_handoff(**options)
    assert _bits(result.membership_pass.memberships) == case["membership_bits"]
    assert _bits(result.membership_pass.objective) == [case["objective_bits"]]
    assert result.registration is not None
    assert result.registration.transform.shape == (6,)
    assert len(result.registration.iterations) <= options["max_iterations"]
    assert result.registration.stop_reason == "iteration_limit"
    assert len(result.registration.iterations) == options["max_iterations"]
    # registerAllShapes resets lambda to INIT_LAMBDA=1 before its first update;
    # this accepted source update grows it by 1.5.
    assert result.registration.iterations[0].lambda_value == np.float32(1.5)
    assert not np.array_equal(result.transform, np.asarray(options["initial_transform"], dtype=np.float32))
    assert result.rotation.shape == (3, 3)
    assert len(result.shape_bounds) == 11
    assert any(item.support_count > 0 for item in result.shape_bounds[1:])
    assert result.membership_pass.memberships.shape == data["initial_memberships"].shape
    output_membership_image = nib.Nifti1Image(result.membership_pass.memberships, subject_affine)
    assert nib.aff2axcodes(output_membership_image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(output_membership_image.affine, subject_affine, rtol=0, atol=0)
    for item in result.shape_bounds:
        assert all(0 <= item.lower_xyz[a] <= data["images"].shape[a] for a in range(3))
        assert all(0 <= item.upper_xyz[a] < data["images"].shape[a] for a in range(3)) or item.support_count == 0

    # Confirm the registration receives the just-computed pass output and that
    # wrapper initialization/reset parameters match the standalone source API.
    initial_bounds = tuple(
        compute_toad_shape_bounds(
            options["atlas_shapes"][..., k],
            image_dimensions_xyz=data["images"].shape[:3],
            image_resolutions_mm=options["image_resolutions_mm"],
            template_resolutions_mm=options["template_resolutions_mm"],
            transform=options["initial_transform"],
        )
        for k in range(11)
    )
    low = np.asarray([item.lower_xyz for item in initial_bounds], dtype=np.int32)
    high = np.asarray([item.upper_xyz for item in initial_bounds], dtype=np.int32)
    direct = register_toads_atlas_rigid(
        expected_membership_from_java,
        options["atlas_shapes"],
        image_resolutions_mm=options["image_resolutions_mm"],
        template_resolutions_mm=options["template_resolutions_mm"],
        initial_transform=options["initial_transform"],
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        registered_classes=np.asarray([False] + [True] * 10),
        max_iterations=options["max_iterations"],
        energy_tolerance=options["energy_tolerance"],
        minimum_lambda=options["minimum_lambda"],
    )
    initial_membership_registration = register_toads_atlas_rigid(
        data["initial_memberships"],
        options["atlas_shapes"],
        image_resolutions_mm=options["image_resolutions_mm"],
        template_resolutions_mm=options["template_resolutions_mm"],
        initial_transform=options["initial_transform"],
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        registered_classes=np.asarray([False] + [True] * 10),
        max_iterations=options["max_iterations"],
        energy_tolerance=options["energy_tolerance"],
        minimum_lambda=options["minimum_lambda"],
    )
    assert abs(float(result.registration.energy - initial_membership_registration.energy)) > 1e-4
    assert np.linalg.norm(
        result.registration.transform - initial_membership_registration.transform
    ) > 1e-4
    np.testing.assert_array_equal(result.registration.transform, direct.transform)
    assert _bits(result.registration.energy) == _bits(direct.energy)
    assert len(result.registration.iterations) == len(direct.iterations)
    for actual, expected in zip(result.registration.iterations, direct.iterations):
        np.testing.assert_array_equal(actual.transform, expected.transform)
        assert _bits(actual.energy) == _bits(expected.energy)
        assert _bits(actual.lambda_value) == _bits(expected.lambda_value)
        assert _bits(actual.relative_change) == _bits(expected.relative_change)
        assert actual.accepted == expected.accepted
        assert actual.rejected_trials == expected.rejected_trials
    expected_final_bounds = tuple(
        compute_toad_shape_bounds(
            options["atlas_shapes"][..., k],
            image_dimensions_xyz=data["images"].shape[:3],
            image_resolutions_mm=options["image_resolutions_mm"],
            template_resolutions_mm=options["template_resolutions_mm"],
            transform=tuple(float(value) for value in result.transform),
        )
        for k in range(11)
    )
    assert result.shape_bounds == expected_final_bounds
    np.testing.assert_allclose(result.rotation @ result.rotation.T, np.eye(3), atol=2e-6)
    for key, value in data.items():
        np.testing.assert_array_equal(value, before[key])
    for key, before_value in statistical_inputs_before.items():
        np.testing.assert_array_equal(options[key], before_value)


def test_align_shapes_false_bypasses_registration_and_preserves_transform() -> None:
    data, parameters, _, _, options, subject_affine = _arguments()
    options["align_shapes"] = False
    result = toads_default_atlas_registration_handoff(**options)
    assert result.registration is None
    np.testing.assert_array_equal(
        result.transform, np.asarray(options["initial_transform"], dtype=np.float32)
    )
    np.testing.assert_allclose(result.rotation @ result.rotation.T, np.eye(3), atol=2e-6)
    assert len(result.shape_bounds) == 11
    assert result.membership_pass.memberships.shape == data["initial_memberships"].shape
    output_membership_image = nib.Nifti1Image(result.membership_pass.memberships, subject_affine)
    assert nib.aff2axcodes(output_membership_image.affine) == ("R", "A", "S")


def test_default_registration_rejects_atlas_grid_mismatch() -> None:
    _, _, _, _, options, _ = _arguments()
    options["atlas_shapes"] = options["atlas_shapes"][..., :-1]
    try:
        toads_default_atlas_registration_handoff(**options)
    except ValueError as error:
        assert "atlas_shapes" in str(error)
    else:
        raise AssertionError("mismatched atlas shapes should reject")
