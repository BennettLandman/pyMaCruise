from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_shapes import (
    accumulate_lesion_registration_coefficients,
    lesion_registration_voxel_contributions,
    lesion_registration_volume_coefficients,
    optimize_lesion_shapes_scaled,
    register_lesion_shapes_scaled_step,
    register_lesion_shapes_pyramid,
    sample_lesion_shape_prior,
)
from pymacruise.toads_registration import trilinear_prior_derivatives


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_shape_sampling"
FIXTURE = FIXTURE_DIR / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionShapeSamplingReferenceHarness.java"
HARNESS_SHA256 = "b5e32e1e77eb66c5edca71b275160963d9d1590d2b77d3e5ccba6a550c78e075"
REGISTRATION_FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_registration_voxel"
REGISTRATION_HARNESS = (
    Path(__file__).parent / "java_reference" / "LesionRegistrationVoxelReferenceHarness.java"
)
REGISTRATION_HARNESS_SHA256 = "68f983a51b6ab834bcc5fcdcf9787170e32a20ac426584d0ae1aeaebc45dbb5e"
SCALED_FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_scaled_shape_sampling"
SCALED_HARNESS = (
    Path(__file__).parent / "java_reference" / "LesionScaledShapeSamplingReferenceHarness.java"
)
SCALED_HARNESS_SHA256 = "379964d519202675bd51fa5a50d561fed22e1e22a9848345f02fb9acde5cc639"
SCALED_FIXTURE_SHA256 = "3c48d9c7cdc57bd652f10b5bdf5c48021ca562877eaaf164b94abb04f49c5b32"
SCALED_REGISTRATION_FIXTURE_DIR = (
    Path(__file__).parent / "data" / "synthetic" / "lesion_scaled_registration_voxel"
)
SCALED_REGISTRATION_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionScaledRegistrationVoxelReferenceHarness.java"
)
SCALED_REGISTRATION_HARNESS_SHA256 = "0fc4d37c5e11fee5a1ef0c491c7a46e1caab5d0e2e64fcdc6e02bbe658de552f"
SCALED_REGISTRATION_FIXTURE_SHA256 = "b00b3abdb5845f829e1cfd22d6590e9c47ea4e10ffa78f3f39745bd97b5a1d2b"
VOLUME_FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_registration_volume"
VOLUME_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionRegistrationVolumeAccumulationReferenceHarness.java"
)
VOLUME_HARNESS_SHA256 = "a2a018ca0c82d575e4bcd195639ffacea7f5b01c6530df87737b4de68d86627c"
VOLUME_FIXTURE_SHA256 = "1320f3f0edc560e8aeac5379708cf0ea392fe6cf3e353b1da6f873b9655f2620"
SCALED_VOLUME_FIXTURE_DIR = (
    Path(__file__).parent / "data" / "synthetic" / "lesion_scaled_registration_volume"
)
SCALED_VOLUME_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionScaledRegistrationVolumeReferenceHarness.java"
)
SCALED_VOLUME_HARNESS_SHA256 = "b624304ff80beafba2500df42abe1e6b5c9309ac4a6a3ef0efa5b28a0eceb187"
SCALED_VOLUME_FIXTURE_SHA256 = "bb20b83650b7062a3c29e1f6c2c3817d01dd90e1a61cd6f0c1257c3aa26a3779"
SCALED_STEP_FIXTURE_DIR = (
    Path(__file__).parent / "data" / "synthetic" / "lesion_scaled_registration_step"
)
SCALED_STEP_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionScaledRegistrationStepReferenceHarness.java"
)
SCALED_STEP_HARNESS_SHA256 = "07b00ea6032b30a1311547123b282de2f84e61ddbf53a67c9a9311f9ab7514b5"
SCALED_STEP_FIXTURE_SHA256 = "8028db4ea1939167b8c572e30eb22b6f2cdac472ab495b188a27bddf48853038"
SCALED_OPTIMIZER_FIXTURE_DIR = (
    Path(__file__).parent / "data" / "synthetic" / "lesion_scaled_registration_optimizer"
)
SCALED_OPTIMIZER_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionScaledRegistrationOptimizerReferenceHarness.java"
)
SCALED_OPTIMIZER_HARNESS_SHA256 = "66fc71104d03f22562bffc94c9e8f947d3edf38d2d18d717467db4e76c761bd0"
SCALED_OPTIMIZER_FIXTURE_SHA256 = "ed9fd9b190a7b6654c24e2f95d7c22052998427b4d51c3bacb0bbd1542ab0794"
LESION_PYRAMID_FIXTURE_DIR = (
    Path(__file__).parent / "data" / "synthetic" / "lesion_registration_pyramid"
)
LESION_PYRAMID_HARNESS = (
    Path(__file__).parent
    / "java_reference"
    / "LesionRegistrationPyramidReferenceHarness.java"
)
LESION_PYRAMID_HARNESS_SHA256 = "3efb5516161b59a6b7c74e3d5e084fb1a623c28982a6af58b231f3d2b5c6a4ec"
LESION_PYRAMID_FIXTURE_SHA256 = "c76c98d37a33eb492601e7674927535da02728565c629f89be29ffaaf900d32f"


def _shape_prior(dimensions: tuple[int, int, int], class_index: int) -> np.ndarray:
    offset = np.float32(0.12 if class_index == 0 else 0.23)
    shape = np.empty(dimensions, dtype=np.float32)
    for x in range(dimensions[0]):
        for y in range(dimensions[1]):
            for z in range(dimensions[2]):
                value = np.float32(offset + np.float32(np.float32(0.017) * x))
                value = np.float32(value + np.float32(np.float32(0.023) * y))
                value = np.float32(value + np.float32(np.float32(0.031) * z))
                value = np.float32(
                    value + np.float32(np.float32(np.float32(0.004) * x) * y)
                )
                shape[x, y, z] = value
    return shape


def test_lesion_shape_sampler_matches_selected_java_formula_harness() -> None:
    reference = json.loads(FIXTURE.read_text())
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for case in reference["cases"]:
        sample = sample_lesion_shape_prior(
            _shape_prior(tuple(reference["template_dims"]), case["class"]),
            tuple(case["point"]),
            image_dimensions_xyz=tuple(reference["image_dims"]),
            image_resolutions_mm=tuple(reference["image_resolutions"]),
            template_resolutions_mm=tuple(reference["template_resolutions"]),
            transform=tuple(case["transform"]),
            class_index=case["class"],
            bbox_low_xyz=tuple(case["bbox_low"]),
            bbox_high_xyz=tuple(case["bbox_high"]),
        )
        assert sample.value == np.float32(case["value"]), case["name"]
        if case["mapped"] is None:
            assert sample.template_xyz is None
        else:
            np.testing.assert_array_equal(
                np.asarray(sample.template_xyz, dtype=np.float32),
                np.asarray(case["mapped"], dtype=np.float32),
                err_msg=case["name"],
            )


def test_lesion_shape_none_mode_directly_samples_ras_loaded_nifti_array(tmp_path: Path) -> None:
    values = np.arange(4 * 5 * 6, dtype=np.float32).reshape((4, 5, 6)) / 100.0
    affine = np.diag([0.8, 1.2, 1.7, 1.0])
    path = tmp_path / "shape.nii.gz"
    nib.save(nib.Nifti1Image(values, affine), path)
    loaded = nib.load(str(path))
    sample = sample_lesion_shape_prior(
        np.asarray(loaded.dataobj, dtype=np.float32),
        (2, 3, 4),
        image_dimensions_xyz=(4, 5, 6),
        image_resolutions_mm=(0.8, 1.2, 1.7),
        template_resolutions_mm=(0.8, 1.2, 1.7),
        transform=(0, 0, 0, 0, 0, 0),
        class_index=2,
        bbox_low_xyz=(0, 0, 0),
        bbox_high_xyz=(4, 5, 6),
        transform_mode="none",
    )
    assert nib.aff2axcodes(loaded.affine) == ("R", "A", "S")
    assert sample.value == values[2, 3, 4]
    assert sample.template_xyz is None


def test_lesion_shape_sampler_rejects_unsupported_modes_and_bad_maps() -> None:
    arguments = {
        "image_dimensions_xyz": (4, 4, 4),
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "transform": (0, 0, 0, 0, 0, 0),
        "class_index": 1,
        "bbox_low_xyz": (0, 0, 0),
        "bbox_high_xyz": (4, 4, 4),
    }
    with pytest.raises(ValueError, match="only 'none' or 'rigid'"):
        sample_lesion_shape_prior(
            np.ones((4, 4, 4), dtype=np.float32),
            (1, 1, 1),
            transform_mode="deformable",
            **arguments,
        )
    with pytest.raises(ValueError, match="match the input grid"):
        sample_lesion_shape_prior(
            np.ones((3, 3, 3), dtype=np.float32),
            (1, 1, 1),
            transform_mode="none",
            **arguments,
        )
    with pytest.raises(ValueError, match="finite"):
        sample_lesion_shape_prior(
            np.full((4, 4, 4), np.nan, dtype=np.float32),
            (1, 1, 1),
            **arguments,
        )


def test_lesion_shape_sampler_exposes_selected_java_prior_gradient() -> None:
    derivative_reference = json.loads(
        (
            Path(__file__).parent
            / "data"
            / "synthetic"
            / "toads_registration_energy"
            / "java_prior_derivative_reference.json"
        ).read_text()
    )
    x, y, z = np.indices((5, 5, 5), dtype=np.int32)
    shape = (
        ((7 * x + 3 * y + 5 * z + 2 * x * y + 4 * y * z) % 23)
        / np.float32(22.0)
    ).astype(np.float32)
    shared = {
        "image_dimensions_xyz": (5, 5, 5),
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "class_index": 1,
        "bbox_low_xyz": (0, 0, 0),
        "bbox_high_xyz": (5, 5, 5),
    }

    fractional = sample_lesion_shape_prior(
        shape,
        (0, 2, 3),
        transform=(0.0, 0.0, 0.0, 0.25, -0.5, -0.25),
        **shared,
    )
    np.testing.assert_array_equal(
        fractional.template_xyz,
        np.asarray((0.25, 1.5, 2.75), dtype=np.float32),
    )
    np.testing.assert_array_equal(
        fractional.template_gradient_xyz,
        np.asarray(derivative_reference["derivatives"][0], dtype=np.float32),
    )

    upper_boundary = sample_lesion_shape_prior(
        shape,
        (3, 3, 3),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        **shared,
    )
    np.testing.assert_array_equal(
        upper_boundary.template_gradient_xyz,
        np.asarray(derivative_reference["derivatives"][1], dtype=np.float32),
    )

    for point, transform, mapped_point, reference_index in (
        ((3, 1, 1), (0, 0, 0, 0.25, 0, 0), (3.25, 1.0, 1.0), 2),
        ((0, 1, 1), (0, 0, 0, -0.1, 0, 0), (-0.1, 1.0, 1.0), 3),
    ):
        outside = sample_lesion_shape_prior(
            shape,
            point,
            transform=transform,
            **shared,
        )
        assert outside.value == np.float32(0.0)
        assert outside.template_gradient_xyz is None
        np.testing.assert_array_equal(
            trilinear_prior_derivatives(shape, mapped_point),
            np.asarray(derivative_reference["derivatives"][reference_index], dtype=np.float32),
        )

    excluded = sample_lesion_shape_prior(
        shape,
        (0, 2, 3),
        image_dimensions_xyz=shared["image_dimensions_xyz"],
        image_resolutions_mm=shared["image_resolutions_mm"],
        template_resolutions_mm=shared["template_resolutions_mm"],
        transform=(0, 0, 0, 0, 0, 0),
        class_index=1,
        bbox_low_xyz=(1, 0, 0),
        bbox_high_xyz=(5, 5, 5),
    )
    assert excluded.template_xyz is None
    assert excluded.template_gradient_xyz is None

    zero_prior = sample_lesion_shape_prior(
        np.zeros((5, 5, 5), dtype=np.float32),
        (2, 2, 2),
        transform=(0, 0, 0, 0, 0, 0),
        **shared,
    )
    assert zero_prior.value == np.float32(0.0)
    assert zero_prior.template_gradient_xyz is None


def test_lesion_registration_voxel_matches_selected_java_formula_fixture() -> None:
    assert hashlib.sha256(REGISTRATION_HARNESS.read_bytes()).hexdigest() == (
        REGISTRATION_HARNESS_SHA256
    )
    reference = json.loads(
        (REGISTRATION_FIXTURE_DIR / "java_reference.json").read_text()
    )

    def affine_map(
        base: float, dx: float, dy: float, dz: float
    ) -> np.ndarray:
        result = np.empty((6, 6, 6), dtype=np.float32)
        for x in range(6):
            for y in range(6):
                for z in range(6):
                    value = np.float32(base)
                    value = np.float32(value + np.float32(np.float32(dx) * x))
                    value = np.float32(value + np.float32(np.float32(dy) * y))
                    value = np.float32(value + np.float32(np.float32(dz) * z))
                    result[x, y, z] = value
        return result

    shapes = (
        np.full((6, 6, 6), 0.9, dtype=np.float32),
        affine_map(0.125, 0.05, 0.025, 0.0125),
        affine_map(0.1, 0.025, 0.05, 0.0125),
        np.zeros((6, 6, 6), dtype=np.float32),
    )
    bounds = tuple((((0, 0, 0)), ((6, 6, 6))) for _ in shapes)
    result = lesion_registration_voxel_contributions(
        shapes,
        (0.9, 0.5, 0.75, 0.6),
        (3, 3, 3),
        image_dimensions_xyz=(6, 6, 6),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        class_bounds_xyz=bounds,
    )

    np.testing.assert_array_equal(
        result.shape_priors, np.asarray(reference["priors"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        result.prior_gradients_xyz,
        np.asarray(reference["prior_gradients_xyz"], dtype=np.float32),
    )
    expected_jacobian = np.asarray(
        ((0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1)),
        dtype=np.float32,
    )
    np.testing.assert_array_equal(result.coordinate_jacobian, expected_jacobian)
    assert result.coefficients.cost == np.float32(reference["cost"])
    assert result.coefficients.norm == np.float32(reference["norm"])
    np.testing.assert_array_equal(
        result.coefficients.gradient,
        np.asarray(reference["gradient"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.coefficients.hessian,
        np.asarray(reference["hessian"], dtype=np.float32),
    )

    excluded_bounds = list(bounds)
    excluded_bounds[1] = ((4, 0, 0), (6, 6, 6))
    excluded = lesion_registration_voxel_contributions(
        shapes,
        (0.9, 0.5, 0.75, 0.6),
        (3, 3, 3),
        image_dimensions_xyz=(6, 6, 6),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        class_bounds_xyz=excluded_bounds,
    )
    assert excluded.shape_priors[1] == np.float32(0.0)
    assert excluded.prior_gradients_xyz[1].tolist() == [0.0, 0.0, 0.0]


def test_lesion_registration_voxel_rejects_class_count_mismatches() -> None:
    shapes = (np.ones((4, 4, 4), dtype=np.float32),) * 2
    with pytest.raises(ValueError, match="one membership and shape prior"):
        lesion_registration_voxel_contributions(
            shapes,
            (0.5,),
            (1, 1, 1),
            image_dimensions_xyz=(4, 4, 4),
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=(0, 0, 0, 0, 0, 0),
            class_bounds_xyz=(((0, 0, 0), (4, 4, 4)),) * 2,
        )
    with pytest.raises(ValueError, match="one bound pair per class"):
        lesion_registration_voxel_contributions(
            shapes,
            (0.5, 0.5),
            (1, 1, 1),
            image_dimensions_xyz=(4, 4, 4),
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=(0, 0, 0, 0, 0, 0),
            class_bounds_xyz=(((0, 0, 0), (4, 4, 4)),),
        )


def test_lesion_shape_sampler_matches_java_registration_pyramid_scales() -> None:
    assert hashlib.sha256(SCALED_HARNESS.read_bytes()).hexdigest() == (
        SCALED_HARNESS_SHA256
    )
    reference_path = SCALED_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == SCALED_FIXTURE_SHA256
    reference = json.loads(reference_path.read_text())
    for case in reference["cases"]:
        dimensions = tuple(case["template_dims"])
        shape = np.empty(dimensions, dtype=np.float32)
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    value = np.float32(0.1)
                    value = np.float32(value + np.float32(np.float32(0.01) * x))
                    value = np.float32(value + np.float32(np.float32(0.02) * y))
                    value = np.float32(value + np.float32(np.float32(0.03) * z))
                    value = np.float32(
                        value + np.float32(np.float32(np.float32(0.001) * x) * y)
                    )
                    shape[x, y, z] = value
        sample = sample_lesion_shape_prior(
            shape,
            tuple(case["point"]),
            image_dimensions_xyz=tuple(case["image_dims"]),
            image_resolutions_mm=tuple(case["image_resolutions"]),
            template_resolutions_mm=tuple(case["template_resolutions"]),
            transform=tuple(case["transform"]),
            class_index=1,
            bbox_low_xyz=(0, 0, 0),
            bbox_high_xyz=tuple(case["image_dims"]),
            scale=case["scale"],
            image_center_xyz=tuple(case["image_center"]),
            template_center_xyz=tuple(case["template_center"]),
        )
        np.testing.assert_array_equal(
            np.asarray(sample.template_xyz, dtype=np.float32),
            np.asarray(case["mapped"], dtype=np.float32),
            err_msg=f"{case['scale']} {case['name']} mapped coordinates",
        )
        assert sample.value == np.float32(case["value"]), (case["scale"], case["name"])
        if case["value"] > 0:
            np.testing.assert_array_equal(
                sample.template_gradient_xyz,
                np.asarray(case["derivative"], dtype=np.float32),
                err_msg=f"{case['scale']} {case['name']} derivatives",
            )
        else:
            assert sample.template_gradient_xyz is None
            np.testing.assert_array_equal(
                trilinear_prior_derivatives(shape, case["mapped"]),
                np.asarray(case["derivative"], dtype=np.float32),
                err_msg=f"{case['scale']} {case['name']} outside derivatives",
            )


def test_lesion_shape_pyramid_sampling_requires_full_resolution_centers() -> None:
    with pytest.raises(ValueError, match="explicit full-resolution centers"):
        sample_lesion_shape_prior(
            np.ones((3, 3, 3), dtype=np.float32),
            (1, 1, 1),
            image_dimensions_xyz=(2, 2, 2),
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=(0, 0, 0, 0, 0, 0),
            class_index=1,
            bbox_low_xyz=(0, 0, 0),
            bbox_high_xyz=(2, 2, 2),
            scale=2,
        )


def test_lesion_scaled_registration_voxel_matches_java_formula_fixture() -> None:
    assert hashlib.sha256(SCALED_REGISTRATION_HARNESS.read_bytes()).hexdigest() == (
        SCALED_REGISTRATION_HARNESS_SHA256
    )
    reference_path = SCALED_REGISTRATION_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        SCALED_REGISTRATION_FIXTURE_SHA256
    )
    reference = json.loads(reference_path.read_text())
    for case in reference["cases"]:
        template_dims = tuple(int(value) for value in case["template_dims"])
        shapes: list[np.ndarray] = []
        for class_index in range(4):
            shape = np.empty(template_dims, dtype=np.float32)
            for x in range(template_dims[0]):
                for y in range(template_dims[1]):
                    for z in range(template_dims[2]):
                        if class_index == 0:
                            value = np.float32(0.9)
                        elif class_index == 1:
                            value = np.float32(0.1)
                            value = np.float32(value + np.float32(np.float32(0.01) * x))
                            value = np.float32(value + np.float32(np.float32(0.02) * y))
                            value = np.float32(value + np.float32(np.float32(0.03) * z))
                            value = np.float32(
                                value + np.float32(np.float32(np.float32(0.001) * x) * y)
                            )
                        elif class_index == 2:
                            value = np.float32(0.15)
                            value = np.float32(value + np.float32(np.float32(0.008) * x))
                            value = np.float32(value + np.float32(np.float32(0.013) * y))
                            value = np.float32(value + np.float32(np.float32(0.021) * z))
                            value = np.float32(
                                value + np.float32(np.float32(np.float32(0.0005) * x) * z)
                            )
                        else:
                            value = np.float32(0.0)
                        shape[x, y, z] = value
            shapes.append(shape)

        input_dims = tuple(int(value) for value in case["image_dims"])
        bounds = tuple(((0, 0, 0), input_dims) for _ in shapes)
        result = lesion_registration_voxel_contributions(
            shapes,
            (0.9, 0.5, 0.75, 0.6),
            tuple(int(value) for value in case["point"]),
            image_dimensions_xyz=input_dims,
            image_resolutions_mm=tuple(case["image_resolutions"]),
            template_resolutions_mm=tuple(case["template_resolutions"]),
            transform=tuple(case["transform"]),
            class_bounds_xyz=bounds,
            scale=case["scale"],
            image_center_xyz=tuple(case["image_center"]),
            template_center_xyz=tuple(case["template_center"]),
        )
        np.testing.assert_array_equal(
            np.asarray(result.template_coordinates_xyz, dtype=np.float32),
            np.asarray(case["mapped"], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.shape_priors, np.asarray(case["priors"], dtype=np.float32)
        )
        np.testing.assert_array_equal(
            result.prior_gradients_xyz,
            np.asarray(case["prior_gradients_xyz"], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.coordinate_jacobian, np.asarray(case["jacobian"], dtype=np.float32)
        )
        assert result.coefficients.cost == np.float32(case["cost"])
        assert result.coefficients.norm == np.float32(case["norm"])
        np.testing.assert_array_equal(
            result.coefficients.gradient,
            np.asarray(case["gradient"], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.coefficients.hessian,
            np.asarray(case["hessian"], dtype=np.float32),
        )


def test_lesion_registration_volume_accumulation_matches_java_loop_fixture() -> None:
    reference_path = VOLUME_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(VOLUME_HARNESS.read_bytes()).hexdigest() == VOLUME_HARNESS_SHA256
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == VOLUME_FIXTURE_SHA256
    reference = json.loads(reference_path.read_text())
    nx, ny, nz = (int(value) for value in reference["dimensions_xyz"])
    classes = int(reference["classes"])
    memberships = np.empty((nx, ny, nz, classes), dtype=np.float32)
    priors = np.empty_like(memberships)
    gradients = np.empty(memberships.shape + (3,), dtype=np.float32)
    jacobians = np.empty((nx, ny, nz, 3, 6), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for k in range(classes):
                    value = np.float32(0.21)
                    value = np.float32(value + np.float32(np.float32(0.09) * k))
                    value = np.float32(value + np.float32(np.float32(0.017) * x))
                    value = np.float32(value + np.float32(np.float32(0.011) * y))
                    memberships[x, y, z, k] = np.float32(
                        value + np.float32(np.float32(0.007) * z)
                    )
                    if k == 0:
                        priors[x, y, z, k] = np.float32(0.88)
                    elif (x + 2 * y + z + k) % 5 == 0:
                        priors[x, y, z, k] = np.float32(0.0)
                    else:
                        value = np.float32(0.11)
                        value = np.float32(value + np.float32(np.float32(0.04) * k))
                        value = np.float32(value + np.float32(np.float32(0.003) * x))
                        value = np.float32(value - np.float32(np.float32(0.002) * y))
                        priors[x, y, z, k] = np.float32(
                            value + np.float32(np.float32(0.001) * z)
                        )
                    gradients[x, y, z, k] = (
                        np.float32(np.float32(0.001) * (k + x + 1)),
                        np.float32(np.float32(0.0013) * (k + y + 2)),
                        np.float32(np.float32(0.0007) * (k + z + 3)),
                    )
                for axis in range(3):
                    for parameter in range(6):
                        value = np.float32(
                            np.float32(np.float32(0.012) * (axis + 1)) * (parameter + 1)
                        )
                        value = np.float32(value + np.float32(np.float32(0.0001) * x))
                        value = np.float32(value - np.float32(np.float32(0.0002) * y))
                        jacobians[x, y, z, axis, parameter] = np.float32(
                            value + np.float32(np.float32(0.0003) * z)
                        )

    results = {}
    for case in reference["cases"]:
        if case["reverse"]:
            continue
        result = accumulate_lesion_registration_coefficients(
            memberships,
            priors,
            gradients,
            jacobians,
            offset=int(case["offset"]),
            subsample=int(case["subsample"]),
        )
        results[case["name"]] = result
        assert result.cost == np.float32(case["cost"]), case["name"]
        assert result.norm == np.float32(case["norm"]), case["name"]
        np.testing.assert_array_equal(result.gradient, np.asarray(case["gradient"], dtype=np.float32))
        np.testing.assert_array_equal(result.hessian, np.asarray(case["hessian"], dtype=np.float32))
        assert result.energy == np.float32(case["energy"]), case["name"]
        np.testing.assert_array_equal(
            result.normalized_gradient,
            np.asarray(case["normalized_gradient"], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.normalized_hessian,
            np.asarray(case["normalized_hessian"], dtype=np.float32),
        )
        assert result.visited_voxels == case["visited_voxels"]
        assert result.contributing_class_terms == case["contributing_class_terms"]

    # The reference includes a reverse-order control to demonstrate that this
    # nonlinear float32 reduction is sensitive to preserving the source order.
    dense = results["dense"]
    reverse = next(case for case in reference["cases"] if case["name"] == "reverse_control")
    assert dense.cost != np.float32(reverse["cost"])
    assert dense.gradient[3] != np.float32(reverse["gradient"][3])


def test_lesion_registration_volume_accumulation_validates_grids_and_source_zero_norm() -> None:
    memberships = np.ones((1, 1, 1, 2), dtype=np.float32)
    priors = np.zeros_like(memberships)
    gradients = np.zeros(memberships.shape + (3,), dtype=np.float32)
    jacobians = np.zeros((1, 1, 1, 3, 6), dtype=np.float32)
    empty = accumulate_lesion_registration_coefficients(
        memberships, priors, gradients, jacobians
    )
    assert empty.cost == empty.norm == np.float32(0.0)
    assert np.isnan(empty.energy)
    assert np.isnan(empty.normalized_gradient).all()
    assert empty.visited_voxels == 1
    assert empty.contributing_class_terms == 0
    with pytest.raises(ValueError, match="positive integer"):
        accumulate_lesion_registration_coefficients(
            memberships, priors, gradients, jacobians, subsample=0
        )
    with pytest.raises(ValueError, match="coordinate_jacobians"):
        accumulate_lesion_registration_coefficients(
            memberships, priors, gradients, np.zeros((1, 1, 1, 3, 5), dtype=np.float32)
        )


def test_scaled_lesion_registration_volume_coefficients_match_java_fixture() -> None:
    reference_path = SCALED_VOLUME_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(SCALED_VOLUME_HARNESS.read_bytes()).hexdigest() == (
        SCALED_VOLUME_HARNESS_SHA256
    )
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        SCALED_VOLUME_FIXTURE_SHA256
    )
    reference = json.loads(reference_path.read_text())
    transform = tuple(reference["transform"])
    classes = int(reference["classes"])
    for case in reference["cases"]:
        input_dims = tuple(int(value) for value in case["input_dims"])
        template_dims = tuple(int(value) for value in case["template_dims"])
        maps = []
        for class_index in range(classes):
            shape = np.empty(template_dims, dtype=np.float32)
            for x in range(template_dims[0]):
                for y in range(template_dims[1]):
                    for z in range(template_dims[2]):
                        if class_index == 0:
                            shape[x, y, z] = np.float32(0.9)
                        else:
                            value = np.float32(0.1)
                            value = np.float32(
                                value + np.float32(np.float32(0.01) * class_index)
                            )
                            value = np.float32(value + np.float32(np.float32(0.013) * x))
                            value = np.float32(value + np.float32(np.float32(0.019) * y))
                            value = np.float32(value + np.float32(np.float32(0.027) * z))
                            shape[x, y, z] = np.float32(
                                value + np.float32(np.float32(np.float32(0.001) * x) * y)
                            )
            maps.append(shape)
        memberships = np.empty(input_dims + (classes,), dtype=np.float32)
        for x in range(input_dims[0]):
            for y in range(input_dims[1]):
                for z in range(input_dims[2]):
                    for class_index in range(classes):
                        value = np.float32(0.23)
                        value = np.float32(
                            value + np.float32(np.float32(0.1) * class_index)
                        )
                        value = np.float32(value + np.float32(np.float32(0.017) * x))
                        value = np.float32(value + np.float32(np.float32(0.011) * y))
                        memberships[x, y, z, class_index] = np.float32(
                            value + np.float32(np.float32(0.007) * z)
                        )
        result = lesion_registration_volume_coefficients(
            maps,
            memberships,
            image_dimensions_xyz=input_dims,
            image_resolutions_mm=tuple(reference["image_resolutions"]),
            template_resolutions_mm=tuple(reference["template_resolutions"]),
            transform=transform,
            scale=case["scale"],
            image_center_xyz=tuple(reference["image_center"]),
            template_center_xyz=tuple(reference["template_center"]),
            offset=0,
            subsample=int(reference["subsample"]),
        )
        assert result.cost == np.float32(case["cost"]), case["scale"]
        assert result.norm == np.float32(case["norm"]), case["scale"]
        np.testing.assert_array_equal(result.gradient, np.asarray(case["gradient"], dtype=np.float32))
        np.testing.assert_array_equal(result.hessian, np.asarray(case["hessian"], dtype=np.float32))
        assert result.energy == np.float32(case["energy"]), case["scale"]
        np.testing.assert_array_equal(
            result.normalized_gradient,
            np.asarray(case["normalized_gradient"], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.normalized_hessian,
            np.asarray(case["normalized_hessian"], dtype=np.float32),
        )
        assert result.visited_voxels == case["visited_voxels"]
        assert result.contributing_class_terms == case["contributing_class_terms"]


def test_scaled_lesion_registration_step_matches_java_trial_trace() -> None:
    reference_path = SCALED_STEP_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(SCALED_STEP_HARNESS.read_bytes()).hexdigest() == (
        SCALED_STEP_HARNESS_SHA256
    )
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        SCALED_STEP_FIXTURE_SHA256
    )
    reference = json.loads(reference_path.read_text())
    expected_final_outcomes = {1: (True, 0), 2: (True, 8)}
    for case in reference["cases"]:
        scale = int(case["scale"])
        dimensions = tuple(int(size) // scale for size in (12, 12, 12))
        shape = np.empty(dimensions, dtype=np.float32)
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    remainder = (
                        15 * x + 4 * y + 10 * z + 5 * x * y + 8 * y * z + 3 * x * z + 3
                    ) % 17
                    shape[x, y, z] = np.float32(
                        np.float32(0.05)
                        + np.float32(np.float32(remainder) / np.float32(20.0))
                    )
        class0 = np.full(dimensions, np.float32(0.9), dtype=np.float32)
        memberships = np.empty(dimensions + (2,), dtype=np.float32)
        a, b, c, d = (23, 15, 8, 14) if scale == 1 else (13, 7, 5, 3)
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    memberships[x, y, z, 0] = np.float32(0.5)
                    remainder = (a * x + b * y + c * z + d) % 23
                    memberships[x, y, z, 1] = np.float32(
                        np.float32(0.04)
                        + np.float32(np.float32(remainder) / np.float32(25.0))
                    )
        step = register_lesion_shapes_scaled_step(
            (class0, shape),
            memberships,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            image_dimensions_xyz=dimensions,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            scale=scale,
            image_center_xyz=(6.0, 6.0, 6.0),
            template_center_xyz=(6.0, 6.0, 6.0),
            offset=0,
            subsample=2,
            initial_lambda=case["initial_lambda"],
        )
        expected_trials = case["trials"]
        np.testing.assert_array_equal(
            step.trial_transforms,
            np.asarray([trial["transform"] for trial in expected_trials], dtype=np.float32),
            err_msg=f"scale {scale} trial transforms",
        )
        np.testing.assert_array_equal(
            step.trial_energies,
            np.asarray([trial["energy"] for trial in expected_trials], dtype=np.float32),
            err_msg=f"scale {scale} trial energies",
        )
        assert step.result.energy == np.float32(case["energy"])
        np.testing.assert_array_equal(
            step.result.transform,
            np.asarray(case["transform"], dtype=np.float32),
        )
        assert step.result.lambda_value == np.float32(case["lambda"])
        assert step.result.relative_change == np.float32(case["relative_change"])
        assert step.result.accepted is case["accepted"]
        assert step.result.rejected_trials == case["rejected_trials"]
        assert (step.result.accepted, step.result.rejected_trials) == expected_final_outcomes[scale]


def test_scaled_lesion_registration_optimizer_loop_matches_java_history() -> None:
    reference_path = SCALED_OPTIMIZER_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(SCALED_OPTIMIZER_HARNESS.read_bytes()).hexdigest() == (
        SCALED_OPTIMIZER_HARNESS_SHA256
    )
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        SCALED_OPTIMIZER_FIXTURE_SHA256
    )
    reference = json.loads(reference_path.read_text())
    dimensions = (8, 8, 8)
    x, y, z = np.indices(dimensions, dtype=np.int32)
    prior_bits = (9 * x + 17 * y + 19 * z + 6 * x * y + 3 * y * z) % 7
    lesion_shape = np.where(prior_bits < 3, np.float32(0.95), np.float32(0.05))
    class0_shape = np.full(dimensions, np.float32(0.9), dtype=np.float32)
    membership_bits = (3 * x + 12 * y + 13 * z + 15 * x * z + 2 * x * y) % 9
    lesion_membership = np.where(
        membership_bits < 4, np.float32(0.9), np.float32(0.1)
    )
    memberships = np.stack(
        (np.full(dimensions, np.float32(0.5)), lesion_membership), axis=-1
    ).astype(np.float32)
    for case in reference["cases"]:
        result = optimize_lesion_shapes_scaled(
            (class0_shape, lesion_shape),
            memberships,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            image_dimensions_xyz=dimensions,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            scale=1.0,
            image_center_xyz=(4.0, 4.0, 4.0),
            template_center_xyz=(4.0, 4.0, 4.0),
            offset=0,
            subsample=2,
            max_iterations=int(case["max_iterations"]),
            min_energy_change=case["min_energy_change"],
            min_lambda=case["min_lambda"],
        )
        assert len(result.iterations) == len(case["iterations"]), case["name"]
        for actual_step, expected_step in zip(
            result.iterations, case["iterations"], strict=True
        ):
            np.testing.assert_array_equal(
                actual_step.trial_transforms,
                np.asarray(
                    [trial["transform"] for trial in expected_step["trials"]],
                    dtype=np.float32,
                ),
            )
            np.testing.assert_array_equal(
                actual_step.trial_energies,
                np.asarray(
                    [trial["energy"] for trial in expected_step["trials"]],
                    dtype=np.float32,
                ),
            )
            actual = actual_step.result
            np.testing.assert_array_equal(
                actual.transform, np.asarray(expected_step["transform"], dtype=np.float32)
            )
            assert actual.energy == np.float32(expected_step["energy"])
            assert actual.lambda_value == np.float32(expected_step["lambda"])
            assert actual.relative_change == np.float32(expected_step["relative_change"])
            assert actual.accepted is expected_step["accepted"]
            assert actual.rejected_trials == expected_step["rejected_trials"]
        np.testing.assert_array_equal(
            result.transform, np.asarray(case["transform"], dtype=np.float32)
        )
        assert result.energy == np.float32(case["energy"])
        assert result.lambda_value == np.float32(case["lambda"])
        assert result.stop_reason == case["stop_reason"]

    outcomes = {case["name"]: case["stop_reason"] for case in reference["cases"]}
    assert outcomes == {
        "source_stops": "lambda_floor",
        "iteration_limit": "iteration_limit",
        "minimum_change": "min_energy_change",
    }


def test_lesion_registration_pyramid_matches_java_levels_and_handoffs() -> None:
    reference_path = LESION_PYRAMID_FIXTURE_DIR / "java_reference.json"
    assert hashlib.sha256(LESION_PYRAMID_HARNESS.read_bytes()).hexdigest() == (
        LESION_PYRAMID_HARNESS_SHA256
    )
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        LESION_PYRAMID_FIXTURE_SHA256
    )
    reference = json.loads(reference_path.read_text())
    dimensions = (33, 33, 33)
    x, y, z = np.indices(dimensions, dtype=np.int32)
    shape_bits = (3 * x + 5 * y + 7 * z + 11 * x * y + 13 * y * z) % 7
    shape1 = np.where(shape_bits < 3, np.float32(0.95), np.float32(0.05))
    shape0 = np.full(dimensions, np.float32(0.9), dtype=np.float32)
    membership_bits = (7 * x + 11 * y + 5 * z + 13 * x * z + 17 * x * y) % 9
    membership1 = np.where(
        membership_bits < 4, np.float32(0.9), np.float32(0.1)
    )
    memberships = np.stack(
        (np.full(dimensions, np.float32(0.5)), membership1), axis=-1
    ).astype(np.float32)
    for case in reference["cases"]:
        result = register_lesion_shapes_pyramid(
            (shape0, shape1),
            memberships,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            levels=int(case["levels"]),
            offset=0,
            subsample=2,
            max_iterations=1,
        )
        assert len(result.levels) == len(case["stages"]), case["levels"]
        for actual_level, expected_level in zip(
            result.levels, case["stages"], strict=True
        ):
            assert actual_level.scale == expected_level["scale"]
            assert list(actual_level.input_dimensions_xyz) == [
                int(size) for size in expected_level["input_dims"]
            ]
            assert list(actual_level.template_dimensions_xyz) == [
                int(size) for size in expected_level["template_dims"]
            ]
            np.testing.assert_array_equal(
                actual_level.initial_transform,
                np.asarray(expected_level["initial_transform"], dtype=np.float32),
            )
            np.testing.assert_array_equal(
                actual_level.optimization.iterations[0].trial_transforms,
                np.asarray(
                    [trial["transform"] for trial in expected_level["trials"]],
                    dtype=np.float32,
                ),
                err_msg=f"level scale {actual_level.scale} trial transforms",
            )
            expected_energies = np.asarray(
                [np.nan if trial["energy"] is None else trial["energy"]
                 for trial in expected_level["trials"]],
                dtype=np.float32,
            )
            np.testing.assert_array_equal(
                actual_level.optimization.iterations[0].trial_energies,
                expected_energies,
                err_msg=f"level scale {actual_level.scale} trial energies",
            )
            iteration = actual_level.optimization.iterations[0].result
            np.testing.assert_array_equal(
                iteration.transform,
                np.asarray(expected_level["transform"], dtype=np.float32),
            )
            assert iteration.energy == np.float32(expected_level["energy"])
            assert iteration.lambda_value == np.float32(expected_level["lambda"])
            assert iteration.relative_change == np.float32(expected_level["relative_change"])
            assert iteration.accepted is expected_level["accepted"]
            assert iteration.rejected_trials == expected_level["rejected_trials"]
            np.testing.assert_array_equal(
                actual_level.optimization.transform,
                np.asarray(expected_level["transform"], dtype=np.float32),
            )
        np.testing.assert_array_equal(
            result.transform, np.asarray(case["final_transform"], dtype=np.float32)
        )
    assert [stage["scale"] for stage in reference["cases"][1]["stages"]] == [2]
    assert [stage["scale"] for stage in reference["cases"][2]["stages"]] == [4, 2]
    assert [stage["scale"] for stage in reference["cases"][3]["stages"]] == [8, 4, 2]
