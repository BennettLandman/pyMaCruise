from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_registration import (
    exact_membership_registration_coefficients,
    exact_membership_registration_energy,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_registration_sampling import (
    build_toads_registration_pyramid,
    register_toads_atlas_rigid,
    register_toads_atlas_rigid_pyramid,
    rigid_registration_coordinate_jacobians,
    sample_toads_registration_fields,
    stream_toads_registration_coefficients,
    stream_toads_registration_energy,
    source_block_mean_subsample,
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)
from pymacruise.toads_shapes import sample_toad_shape_prior


PYRAMID_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_registration_pyramid"
    / "java_reference.json"
)
PYRAMID_HARNESS_SHA256 = "3893d0fed7f3635ba1002a07c08bcf021475cd72fc95ddde68e8d78da3dc9f32"
RAS_MAPPING_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_preprocessing"
    / "ras_mapping_java_reference.json"
)
RAS_MAPPING_HARNESS_SHA256 = "0c79963fd337053f7a977f44e22202dc89d9738d515821f667611f4752e05128"
TOADS_REGISTRATION_PYRAMID_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_registration_pyramid_optimizer"
)
NONFINITE_TRIAL_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_registration_nonfinite_trial"
    / "java_reference.json"
)


def test_source_block_mean_subsampling_matches_java_harness() -> None:
    reference = json.loads(PYRAMID_FIXTURE.read_text())
    x, y, z, k = np.indices((9, 8, 7, 3), dtype=np.int32)
    values = (((11 * x + 5 * y + 3 * z + 7 * k) % 29) / np.float32(28.0)).astype(np.float32)
    for factor in (2, 4):
        expected = reference[str(factor)]
        actual = source_block_mean_subsample(values, factor)
        assert list(actual.shape) == expected["shape"]
        np.testing.assert_array_equal(actual.reshape(-1), np.asarray(expected["values"], np.float32))


def test_pyramid_subsampling_harness_hash_and_coarse_to_fine_schedule() -> None:
    harness = Path(__file__).parent / "java_reference" / "ToadsPyramidSubsampleHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == PYRAMID_HARNESS_SHA256
    assert [level.scale for level in build_toads_registration_pyramid(
        np.ones((32, 32, 32, 2), np.float32),
        np.ones((32, 32, 32, 2), np.float32),
        levels=1,
    )] == []
    for levels in range(2, 6):
        pyramid = build_toads_registration_pyramid(
            np.ones((32, 32, 32, 2), np.float32),
            np.ones((32, 32, 32, 2), np.float32),
            levels=levels,
        )
        # The source halves integer power-of-two scales at each step.
        expected_scales = [2**power for power in range(levels - 1, 0, -1)]
        assert [level.scale for level in pyramid] == expected_scales
        assert all(level.memberships.shape[:3] == tuple(size // level.scale for size in (32, 32, 32)) for level in pyramid)
        assert all(level.atlas_shapes.shape == level.memberships.shape for level in pyramid)


@pytest.mark.parametrize(
    ("image", "factor", "message"),
    [
        (np.ones((3, 4), np.float32), 2, "3D or 4D"),
        (np.ones((3, 4, 5), np.float32), 0, "positive integer"),
        (np.ones((3, 4, 5), np.float32), 4, "too large"),
        (np.full((3, 4, 5), np.nan, np.float32), 2, "finite"),
    ],
)
def test_source_block_mean_subsampling_rejects_invalid_inputs(
    image: np.ndarray, factor: int, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        source_block_mean_subsample(image, factor)


def _atlas() -> np.ndarray:
    x, y, z, k = np.indices((6, 5, 4, 3), dtype=np.float32)
    return (
        0.11
        + 0.013 * x
        + 0.019 * y
        + 0.027 * z
        + 0.031 * k
        + 0.004 * x * y
        - 0.003 * y * z
    ).astype(np.float32)


def test_transform_dependent_registration_fields_compose_source_kernels() -> None:
    atlas = _atlas()
    dimensions = (4, 3, 3)
    image_resolutions = (1.4, 1.1, 1.7)
    template_resolutions = (0.9, 1.2, 1.5)
    transform = (0.08, -0.04, 0.03, 0.5, -0.25, 0.2)
    lows = np.asarray(((0, 0, 0), (1, 0, 0), (0, 1, 1)), dtype=np.int64)
    highs = np.asarray((dimensions, (4, 3, 3), (3, 3, 3)), dtype=np.int64)

    fields = sample_toads_registration_fields(
        atlas,
        input_dimensions_xyz=dimensions,
        image_resolutions_mm=image_resolutions,
        template_resolutions_mm=template_resolutions,
        transform=transform,
        class_bounds_low_xyz=lows,
        class_bounds_high_xyz=highs,
    )

    assert fields.priors.shape == dimensions + (3,)
    assert fields.prior_gradients_xyz.shape == dimensions + (3, 3)
    assert fields.template_coordinates_xyz.shape == dimensions + (3,)
    for point in ((0, 0, 0), (1, 1, 1), (3, 2, 2)):
        x, y, z = point
        for class_index in range(3):
            expected = sample_toad_shape_prior(
                atlas[..., class_index],
                point,
                image_dimensions_xyz=dimensions,
                image_resolutions_mm=image_resolutions,
                template_resolutions_mm=template_resolutions,
                transform=transform,
                class_index=class_index,
                bbox_low_xyz=tuple(lows[class_index]),
                bbox_high_xyz=tuple(highs[class_index]),
                shape_slope=1.0,
            )
            assert fields.priors[x, y, z, class_index] == expected.value
            if expected.template_xyz is None:
                np.testing.assert_array_equal(
                    fields.prior_gradients_xyz[x, y, z, class_index], np.zeros(3, np.float32)
                )
            else:
                np.testing.assert_array_equal(
                    fields.template_coordinates_xyz[x, y, z], expected.template_xyz
                )
                # Outside derivative support is deliberately zero even where
                # closest-edge value interpolation still returns a value.
                np.testing.assert_array_equal(
                    fields.prior_gradients_xyz[x, y, z, class_index],
                    trilinear_prior_derivatives(atlas[..., class_index], expected.template_xyz),
                )

    # Registration keeps unsloped raw priors. The segmentation helper's
    # default slope recovery is intentionally a different contract.
    point = (1, 1, 1)
    segmentation_sample = sample_toad_shape_prior(
        atlas[..., 1], point,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=image_resolutions,
        template_resolutions_mm=template_resolutions,
        transform=transform,
        class_index=1,
        bbox_low_xyz=tuple(lows[1]),
        bbox_high_xyz=tuple(highs[1]),
    )
    assert fields.priors[point + (1,)] != segmentation_sample.value


def test_ras_signed_spacing_landmarks_match_selected_java_mapping() -> None:
    reference_root = Path(__file__).parent / "java_reference"
    harness = reference_root / "ToadsOrientationMappingHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == RAS_MAPPING_HARNESS_SHA256
    reference = json.loads(RAS_MAPPING_FIXTURE.read_text())
    dims = (7, 6, 5)
    atlas_dims = (11, 10, 9)
    fields = sample_toads_registration_fields(
        np.ones(atlas_dims + (2,), dtype=np.float32),
        input_dimensions_xyz=dims,
        image_resolutions_mm=tuple(reference["signed_resolutions"]),
        template_resolutions_mm=(0.8, 0.9, 1.1),
        transform=(0.0,) * 6,
        class_bounds_low_xyz=np.zeros((2, 3), dtype=np.int64),
        class_bounds_high_xyz=np.tile(np.asarray(dims), (2, 1)),
    )
    for point, expected in zip(((0, 0, 0), (3, 2, 1), (6, 5, 4)), reference["points"], strict=True):
        np.testing.assert_array_equal(
            fields.template_coordinates_xyz[point], np.asarray(expected, dtype=np.float32)
        )


def test_rigid_registration_jacobian_field_matches_single_voxel_kernel() -> None:
    dimensions = (3, 2, 4)
    transform = (0.08, -0.04, 0.03, 0.5, -0.25, 0.2)
    image_res = (1.4, 1.1, 1.7)
    template_res = (0.9, 1.2, 1.5)
    field = rigid_registration_coordinate_jacobians(
        dimensions, transform, image_res, template_res
    )
    center = tuple(np.float32(size / 2.0) for size in dimensions)
    for point in ((0, 0, 0), (1, 1, 2), (2, 1, 3)):
        np.testing.assert_array_equal(
            field[point],
            rigid_image_to_template_jacobian(
                point, transform, center, image_res, template_res
            ),
        )


def test_transform_dependent_registration_fields_reject_malformed_shapes_and_bounds() -> None:
    atlas = np.ones((3, 3, 3, 2), dtype=np.float32)
    common = dict(
        input_dimensions_xyz=(2, 2, 2),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        class_bounds_low_xyz=np.zeros((2, 3), dtype=np.int64),
        class_bounds_high_xyz=np.full((2, 3), 2, dtype=np.int64),
    )
    sample_toads_registration_fields(atlas, **common)

    with pytest.raises(ValueError, match="4D"):
        sample_toads_registration_fields(atlas[..., 0], **common)
    with pytest.raises(ValueError, match=r"shape \(class, 3\)"):
        sample_toads_registration_fields(
            atlas, **{**common, "class_bounds_high_xyz": np.full((3, 2), 2)}
        )


def test_rigid_registration_optimizer_runs_repeatably_on_seeded_synthetic_atlas() -> None:
    atlas_shape = (5, 5, 5, 3)
    x, y, z, class_index = np.indices(atlas_shape, dtype=np.float32)
    atlas = (
        0.18 + 0.035 * x + 0.025 * y + 0.02 * z + 0.03 * class_index + 0.004 * x * y
    ).astype(np.float32)
    dimensions = (4, 4, 4)
    ix, iy, iz = np.indices(dimensions, dtype=np.float32)
    memberships = np.empty(dimensions + (3,), dtype=np.float32)
    memberships[..., 0] = 0.1
    memberships[..., 1] = 0.25 + 0.04 * ix + 0.02 * iy
    memberships[..., 2] = 0.2 + 0.015 * iz + 0.01 * ix
    bounds_low = np.zeros((3, 3), dtype=np.int64)
    bounds_high = np.full((3, 3), 4, dtype=np.int64)
    kwargs = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.0, 0.0, 0.0, 0.05, -0.03, 0.02),
        class_bounds_low_xyz=bounds_low,
        class_bounds_high_xyz=bounds_high,
        max_iterations=2,
    )

    first = register_toads_atlas_rigid(memberships, atlas, **kwargs)
    second = register_toads_atlas_rigid(memberships, atlas, **kwargs)

    assert first.stop_reason == "iteration_limit"
    assert len(first.iterations) == 2
    assert all(item.accepted for item in first.iterations)
    assert all(item.rejected_trials < 11 for item in first.iterations)
    assert np.isfinite(first.energy)
    np.testing.assert_array_equal(first.transform, second.transform)
    assert first.energy == second.energy
    for left, right in zip(first.iterations, second.iterations, strict=True):
        np.testing.assert_array_equal(left.transform, right.transform)
        assert left.energy == right.energy
        assert left.lambda_value == right.lambda_value
        assert left.relative_change == right.relative_change


def test_rigid_registration_optimizer_rejects_empty_positive_prior_support() -> None:
    memberships = np.ones((2, 2, 2, 2), dtype=np.float32)
    atlas = np.zeros((3, 3, 3, 2), dtype=np.float32)
    bounds_low = np.zeros((2, 3), dtype=np.int64)
    bounds_high = np.full((2, 3), 2, dtype=np.int64)
    with pytest.raises(ValueError, match="no positive prior support"):
        register_toads_atlas_rigid(
            memberships,
            atlas,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            initial_transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            class_bounds_low_xyz=bounds_low,
            class_bounds_high_xyz=bounds_high,
        )


def test_rigid_registration_optimizer_honors_outer_stopping_gates() -> None:
    atlas_shape = (4, 4, 4, 2)
    x, y, z, _ = np.indices(atlas_shape, dtype=np.float32)
    atlas = (0.2 + 0.03 * x + 0.02 * y + 0.01 * z).astype(np.float32)
    memberships = np.empty((3, 3, 3, 2), dtype=np.float32)
    mx, my, mz = np.indices((3, 3, 3), dtype=np.float32)
    memberships[..., 0] = 0.1
    memberships[..., 1] = 0.3 + 0.02 * mx + 0.01 * my + 0.01 * mz
    kwargs = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.0, 0.0, 0.0, 0.02, -0.01, 0.01),
        class_bounds_low_xyz=np.zeros((2, 3), dtype=np.int64),
        class_bounds_high_xyz=np.full((2, 3), 3, dtype=np.int64),
        max_iterations=5,
        energy_tolerance=1.0,
    )
    by_energy = register_toads_atlas_rigid(memberships, atlas, **kwargs)
    assert by_energy.stop_reason == "energy_tolerance"
    assert len(by_energy.iterations) == 1

    by_lambda = register_toads_atlas_rigid(
        memberships, atlas, **{**kwargs, "energy_tolerance": 0.0, "minimum_lambda": 2.0}
    )
    assert by_lambda.stop_reason == "lambda_floor"
    assert len(by_lambda.iterations) == 1


def test_streaming_registration_matches_dense_fields_and_coefficients_exactly() -> None:
    atlas = _atlas()
    dimensions = (4, 3, 3)
    image_res = (1.4, 1.1, 1.7)
    template_res = (0.9, 1.2, 1.5)
    transform = (0.08, -0.04, 0.03, 0.5, -0.25, 0.2)
    low = np.asarray(((0, 0, 0), (1, 0, 0), (0, 1, 1)), dtype=np.int64)
    high = np.asarray((dimensions, (4, 3, 3), (3, 3, 3)), dtype=np.int64)
    x, y, z, k = np.indices(dimensions + (3,), dtype=np.int32)
    memberships = (((3 * x + 5 * y + 7 * z + 2 * k) % 17) / np.float32(16.0)).astype(np.float32)
    registered = (True, True, False)
    fields = sample_toads_registration_fields(
        atlas,
        input_dimensions_xyz=dimensions,
        image_resolutions_mm=image_res,
        template_resolutions_mm=template_res,
        transform=transform,
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
    )
    jacobians = rigid_registration_coordinate_jacobians(
        dimensions, transform, image_res, template_res
    )

    for offset, subsample in ((0, 1), (1, 2)):
        dense_energy = exact_membership_registration_energy(
            memberships,
            fields.priors,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        streaming_energy = stream_toads_registration_energy(
            memberships,
            atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        assert streaming_energy == dense_energy

        dense = exact_membership_registration_coefficients(
            memberships,
            fields.priors,
            fields.prior_gradients_xyz,
            jacobians,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        streaming = stream_toads_registration_coefficients(
            memberships,
            atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        for name in ("cost", "norm"):
            assert getattr(streaming, name) == getattr(dense, name)
        for name in ("gradient", "hessian", "norm_gradient", "norm_hessian"):
            np.testing.assert_array_equal(getattr(streaming, name), getattr(dense, name))

        vectorized_energy = vectorized_toads_registration_energy(
            memberships,
            atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        assert vectorized_energy == dense_energy
        vectorized = vectorized_toads_registration_coefficients(
            memberships,
            atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        for name in ("cost", "norm"):
            assert getattr(vectorized, name) == getattr(dense, name)
        for name in ("gradient", "hessian", "norm_gradient", "norm_hessian"):
            np.testing.assert_array_equal(getattr(vectorized, name), getattr(dense, name))


def test_vectorized_registration_matches_scalar_with_large_rigid_offset() -> None:
    atlas = _atlas()
    dimensions = (3, 4, 2)
    image_res = (1.7, 0.8, 2.1)
    template_res = (0.9, 1.2, 0.7)
    transform = (0.21, -0.11, 0.07, 2.0, -1.5, 1.1)
    low = np.asarray(((0, 0, 0), (1, 1, 0), (0, 0, 1)), dtype=np.int64)
    high = np.asarray((dimensions, (3, 4, 1), (3, 4, 2)), dtype=np.int64)
    x, y, z, k = np.indices(dimensions + (3,), dtype=np.int32)
    memberships = (((5 * x + 3 * y + 7 * z + 2 * k) % 13) / np.float32(12.0)).astype(np.float32)
    fields = sample_toads_registration_fields(
        atlas,
        input_dimensions_xyz=dimensions,
        image_resolutions_mm=image_res,
        template_resolutions_mm=template_res,
        transform=transform,
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
    )
    jacobians = rigid_registration_coordinate_jacobians(
        dimensions, transform, image_res, template_res
    )
    for offset, subsample in ((0, 1), (1, 2)):
        dense = exact_membership_registration_coefficients(
            memberships, fields.priors, fields.prior_gradients_xyz, jacobians,
            offset=offset, subsample=subsample,
        )
        scalar = stream_toads_registration_coefficients(
            memberships, atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            offset=offset,
            subsample=subsample,
        )
        vectorized = vectorized_toads_registration_coefficients(
            memberships, atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            offset=offset,
            subsample=subsample,
        )
        for field in ("cost", "norm", "gradient", "hessian", "norm_gradient", "norm_hessian"):
            np.testing.assert_array_equal(getattr(scalar, field), getattr(dense, field))
            np.testing.assert_array_equal(getattr(vectorized, field), getattr(dense, field))


def test_pyramid_registration_runs_scales_and_hands_transform_between_levels() -> None:
    membership_dims = (17, 15, 13)
    atlas_dims = (19, 17, 15)
    fixture_dir = Path(__file__).parent / "data" / "synthetic" / "toads_registration_pyramid_optimizer"
    membership_image = nib.load(fixture_dir / "memberships.nii.gz")
    atlas_image = nib.load(fixture_dir / "atlas_shapes.nii.gz")
    assert nib.aff2axcodes(membership_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(atlas_image.affine) == ("R", "A", "S")
    memberships = np.asanyarray(membership_image.dataobj, dtype=np.float32)
    atlas = np.asanyarray(atlas_image.dataobj, dtype=np.float32)
    assert memberships.shape == membership_dims + (3,)
    assert atlas.shape == atlas_dims + (3,)
    rng = np.random.default_rng(2)
    np.testing.assert_array_equal(
        memberships,
        rng.uniform(0.05, 0.9, size=membership_dims + (3,)).astype(np.float32),
    )
    np.testing.assert_array_equal(
        atlas,
        rng.uniform(0.05, 0.9, size=atlas_dims + (3,)).astype(np.float32),
    )
    initial = np.asarray((0.04, -0.025, 0.015, 0.03, -0.02, 0.01), dtype=np.float32)
    kwargs = dict(
        levels=3,
        image_resolutions_mm=(0.8, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.1, 1.3),
        initial_transform=tuple(float(v) for v in initial),
        max_iterations=2,
    )

    result = register_toads_atlas_rigid_pyramid(memberships, atlas, **kwargs)
    repeated = register_toads_atlas_rigid_pyramid(memberships, atlas, **kwargs)

    assert [level.scale for level in result.levels] == [4, 2]
    assert any(
        iteration.accepted
        for level in result.levels
        for iteration in level.result.iterations
    )
    np.testing.assert_array_equal(result.levels[0].initial_transform, initial)
    np.testing.assert_array_equal(
        result.levels[1].initial_transform, result.levels[0].result.transform
    )
    np.testing.assert_array_equal(result.transform, result.levels[-1].result.transform)
    np.testing.assert_array_equal(result.transform, repeated.transform)
    for left, right in zip(result.levels, repeated.levels, strict=True):
        np.testing.assert_array_equal(left.initial_transform, right.initial_transform)
        np.testing.assert_array_equal(left.result.transform, right.result.transform)
        np.testing.assert_equal(left.result.energy, right.result.energy)
        assert len(left.result.iterations) <= 2
    assert np.isfinite(result.levels[-1].result.energy)


def test_pyramid_registration_matches_selected_java_reference() -> None:
    fixture_dir = TOADS_REGISTRATION_PYRAMID_FIXTURE
    reference = json.loads((fixture_dir / "java_reference.json").read_text())
    memberships = np.asanyarray(
        nib.load(fixture_dir / "memberships.nii.gz").dataobj, dtype=np.float32
    )
    atlas = np.asanyarray(
        nib.load(fixture_dir / "atlas_shapes.nii.gz").dataobj, dtype=np.float32
    )
    initial = tuple(reference["initial_transform"])
    result = register_toads_atlas_rigid_pyramid(
        memberships,
        atlas,
        levels=3,
        image_resolutions_mm=(0.8, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.1, 1.3),
        initial_transform=initial,
        max_iterations=2,
        # This selected-helper reference exercises its independently
        # transcribed all-voxel/MultipleToadAtlas coefficient convention.
        # The actual VandyCRUISE method uses offset-zero step-two sampling and
        # the MultipleToadDeformableAtlas gradient path tested below.
        subsample_factor=1,
        source_faithful=True,
    )

    assert [level.scale for level in result.levels] == [4, 2]
    assert len(result.levels) == len(reference["stages"])
    np.testing.assert_array_equal(result.levels[0].initial_transform, np.asarray(initial, np.float32))
    for actual, expected in zip(result.levels, reference["stages"], strict=True):
        assert actual.scale == expected["scale"]
        assert list(actual.result.transform.shape) == [6]
        np.testing.assert_array_equal(actual.initial_transform, np.asarray(expected["initial_transform"], np.float32))
        np.testing.assert_array_equal(actual.result.transform, np.asarray(expected["transform"], np.float32))
        assert actual.result.energy == np.float32(expected["energy"])
        assert actual.result.iterations[-1].lambda_value == np.float32(expected["lambda"])
        assert actual.result.iterations[-1].relative_change == np.float32(expected["relative_change"])
        assert len(actual.result.iterations) == len(expected["iterations"])
        stop_reason = (
            "energy_tolerance"
            if actual.result.iterations[-1].relative_change <= np.float32(1e-6)
            else "lambda_floor"
            if actual.result.iterations[-1].lambda_value <= np.float32(1e-3)
            else "iteration_limit"
        )
        assert actual.result.stop_reason == stop_reason
        for actual_iteration, expected_iteration in zip(
            actual.result.iterations, expected["iterations"], strict=True
        ):
            np.testing.assert_array_equal(
                actual_iteration.transform,
                np.asarray(expected_iteration["transform"], np.float32),
            )
            assert actual_iteration.energy == np.float32(expected_iteration["energy"])
            assert actual_iteration.lambda_value == np.float32(expected_iteration["lambda"])
            assert actual_iteration.relative_change == np.float32(expected_iteration["relative_change"])
            assert actual_iteration.accepted is expected_iteration["accepted"]
            assert actual_iteration.rejected_trials == expected_iteration["rejected_trials"]
            assert len(expected_iteration["trials"]) == (
                expected_iteration["rejected_trials"]
                + int(expected_iteration["accepted"])
            )


def test_pyramid_registration_matches_actual_vandycruise_entry_point(tmp_path: Path) -> None:
    fixture_dir = TOADS_REGISTRATION_PYRAMID_FIXTURE
    memberships_image = nib.load(fixture_dir / "memberships.nii.gz")
    atlas_image = nib.load(fixture_dir / "atlas_shapes.nii.gz")
    memberships = np.asanyarray(memberships_image.dataobj, dtype=np.float32)
    atlas = np.asanyarray(atlas_image.dataobj, dtype=np.float32)
    assert nib.aff2axcodes(memberships_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(atlas_image.affine) == ("R", "A", "S")
    mem_path = tmp_path / "memberships-f32le.bin"
    atlas_path = tmp_path / "atlas-f32le.bin"
    memberships.tofile(mem_path)
    atlas.tofile(atlas_path)
    image_spacing = memberships_image.header.get_zooms()[:3]
    atlas_spacing = atlas_image.header.get_zooms()[:3]
    initial = (0.04, -0.025, 0.015, 0.03, -0.02, 0.01)
    command = [
        str(Path(__file__).parent / "java_reference" / "run_toads_pyramid_actual_registration_reference.sh"),
        str(mem_path), str(atlas_path),
        *(str(value) for value in memberships.shape[:3]),
        *(str(value) for value in atlas.shape[:3]),
        str(memberships.shape[3]),
        *(str(float(value)) for value in image_spacing),
        *(str(float(value)) for value in atlas_spacing),
        "3", "2", *(str(value) for value in initial),
    ]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    java_result = json.loads(next(line for line in completed.stdout.splitlines() if line.startswith("{")))
    assert java_result["source"].endswith("registerShapesExactPyramid()")

    python_result = register_toads_atlas_rigid_pyramid(
        memberships,
        atlas,
        levels=3,
        image_resolutions_mm=tuple(float(value) for value in image_spacing),
        template_resolutions_mm=tuple(float(value) for value in atlas_spacing),
        initial_transform=initial,
        registered_classes=np.asarray((False, True, True)),
        max_iterations=2,
        allow_nonfinite_trials=True,
    )
    np.testing.assert_array_equal(
        np.asarray(java_result["transform_bits"], dtype=np.int32),
        python_result.transform.view(np.int32),
        err_msg=("actual VandyCRUISE pyramid source differs from the Python port; source trace:\n"
                 + completed.stdout),
    )
    assert np.int32(java_result["lambda_bits"]) == python_result.levels[-1].result.iterations[-1].lambda_value.view(np.int32)
    assert python_result.levels[-1].result.energy.view(np.int32) == np.float32(java_result["final_energy"]).view(np.int32)
    coarse = python_result.levels[0].result
    assert len(coarse.iterations) == 1
    assert coarse.iterations[0].rejected_trials == 11
    assert not coarse.iterations[0].accepted
    assert coarse.energy.view(np.int32) == np.float32(java_result["coarse_nan_trial_energy"]).view(np.int32)
    coarse_log = completed.stdout.split("level 2", maxsplit=1)[0]
    source_initial_match = re.search(r"a: 1\.0\(([^)]+)\)->", coarse_log)
    assert source_initial_match is not None
    source_trial_energies = re.findall(r"\)->\(([^)]+)\)", coarse_log)
    assert len(source_trial_energies) == coarse.iterations[0].rejected_trials
    assert coarse.energy.view(np.int32) == np.float32(source_trial_energies[-1]).view(np.int32)
    coarse_delta = np.float32(coarse.energy - np.float32(source_initial_match.group(1)))
    source_relative_change = np.float32(coarse_delta / coarse.energy)
    assert coarse.iterations[0].relative_change.view(np.int32) == source_relative_change.view(np.int32)
    assert coarse_log.count("trial transform:") == coarse.iterations[0].rejected_trials
    np.testing.assert_array_equal(coarse.transform, np.asarray(initial, np.float32))
    assert " (E:" not in coarse_log
    fine_log = completed.stdout.split("level 2", maxsplit=1)[1]
    assert fine_log.count("trial transform:") == sum(
        iteration.rejected_trials + int(iteration.accepted)
        for iteration in python_result.levels[-1].result.iterations
    )
    assert fine_log.count(" (E:") == sum(
        iteration.accepted for iteration in python_result.levels[-1].result.iterations
    )


def test_pyramid_registration_levels_one_is_source_noop() -> None:
    memberships = np.ones((5, 4, 3, 2), dtype=np.float32)
    atlas = np.ones((7, 6, 5, 2), dtype=np.float32)
    initial = np.asarray((0.0, 0.0, 0.0, 0.2, -0.1, 0.05), dtype=np.float32)
    result = register_toads_atlas_rigid_pyramid(
        memberships,
        atlas,
        levels=1,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=tuple(float(v) for v in initial),
    )
    assert result.levels == ()
    np.testing.assert_array_equal(result.transform, initial)


def test_nonfinite_registration_trial_matches_source_rejection_and_stops_scale() -> None:
    reference = json.loads(NONFINITE_TRIAL_FIXTURE.read_text())
    harness = (
        Path(__file__).parent
        / "java_reference"
        / "ToadsNonfiniteRegistrationTrialReferenceHarness.java"
    )
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == reference["harness_sha256"]
    dimensions = (3, 3, 3)
    memberships = np.full(dimensions + (2,), np.float32(0.5), dtype=np.float32)
    shapes = np.zeros(dimensions + (2,), dtype=np.float32)
    shapes[..., 1] = np.float32(0.5)
    initial = np.zeros(6, dtype=np.float32)
    bounds_low = np.zeros((2, 3), dtype=np.int64)
    bounds_high = np.broadcast_to(np.asarray(dimensions), (2, 3)).copy()

    result = register_toads_atlas_rigid(
        memberships,
        shapes,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=tuple(float(value) for value in initial),
        class_bounds_low_xyz=bounds_low,
        class_bounds_high_xyz=bounds_high,
        registered_classes=(False, True),
        max_iterations=5,
    )

    assert result.stop_reason == "nonfinite_relative_change"
    assert len(result.iterations) == 1
    actual = {
        "accepted": result.iterations[0].accepted,
        "rejected_trials": result.iterations[0].rejected_trials,
        "energy_nan": bool(np.isnan(result.iterations[0].energy)),
        "relative_change_nan": bool(np.isnan(result.iterations[0].relative_change)),
        "outer_continues": bool(result.iterations[0].relative_change > np.float32(1e-6)),
        "transform_unchanged": bool(np.array_equal(result.transform, initial)),
    }
    for key, expected in reference.items():
        if key in actual:
            assert actual[key] == expected, key


def test_scaled_vectorized_registration_matches_scalar_reference_kernels() -> None:
    dimensions = (3, 3, 2)
    template_dimensions = (4, 3, 3)
    image_center = (5.5, 5.0, 3.5)
    template_center = (7.5, 5.5, 4.5)
    image_res = (0.8, 1.2, 1.5)
    template_res = (0.9, 1.1, 1.3)
    scale = 2.0
    transform = (0.14, 0.07, -0.05, 0.2, -0.1, 0.3)
    ax, ay, az = np.indices(template_dimensions, dtype=np.float32)
    atlas = np.empty(template_dimensions + (3,), dtype=np.float32)
    atlas[..., 0] = 0.0
    atlas[..., 1] = 0.2 + 0.04 * ax + 0.03 * ay + 0.02 * az + 0.005 * ax * ay
    atlas[..., 2] = 0.15 + 0.025 * ax + 0.05 * ay + 0.03 * az + 0.004 * ay * az
    mx, my, mz, mk = np.indices(dimensions + (3,), dtype=np.float32)
    memberships = (0.12 + 0.03 * mx + 0.02 * my + 0.01 * mz + 0.015 * mk).astype(np.float32)
    low = np.zeros((3, 3), dtype=np.int64)
    high = np.broadcast_to(np.asarray(dimensions, dtype=np.int64), low.shape).copy()
    priors = np.zeros(dimensions + (3,), dtype=np.float32)
    gradients = np.zeros(dimensions + (3, 3), dtype=np.float32)
    jacobians = np.empty(dimensions + (3, 6), dtype=np.float32)
    for x in range(dimensions[0]):
        for y in range(dimensions[1]):
            for z in range(dimensions[2]):
                jacobians[x, y, z] = rigid_image_to_template_jacobian(
                    (x, y, z), transform, image_center, image_res, template_res,
                    scale=scale,
                )
                for class_index in range(1, 3):
                    sample = sample_toad_shape_prior(
                        atlas[..., class_index],
                        (x, y, z),
                        image_dimensions_xyz=dimensions,
                        image_resolutions_mm=image_res,
                        template_resolutions_mm=template_res,
                        transform=transform,
                        class_index=class_index,
                        bbox_low_xyz=(0, 0, 0),
                        bbox_high_xyz=dimensions,
                        shape_slope=1.0,
                        scale=scale,
                        image_center_xyz=image_center,
                        template_center_xyz=template_center,
                    )
                    assert sample.template_xyz is not None
                    priors[x, y, z, class_index] = sample.value
                    gradients[x, y, z, class_index] = trilinear_prior_derivatives(
                        atlas[..., class_index], sample.template_xyz
                    )

    reference = exact_membership_registration_coefficients(
        memberships, priors, gradients, jacobians
    )
    actual = vectorized_toads_registration_coefficients(
        memberships,
        atlas,
        image_resolutions_mm=image_res,
        template_resolutions_mm=template_res,
        transform=transform,
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        scale=scale,
        image_center_xyz=image_center,
        template_center_xyz=template_center,
    )
    max_absolute_difference = 0.0
    max_difference_field = ""
    for name in ("cost", "norm", "gradient", "hessian", "norm_gradient", "norm_hessian"):
        field_difference = float(
            np.max(np.abs(getattr(actual, name) - getattr(reference, name)))
        )
        if field_difference > max_absolute_difference:
            max_absolute_difference = field_difference
            max_difference_field = name
        np.testing.assert_allclose(
            getattr(actual, name),
            getattr(reference, name),
            rtol=2e-7,
            atol=3e-10,
            err_msg=name,
        )
    assert max_absolute_difference <= 1.5e-8, (max_difference_field, max_absolute_difference)
    reference_energy = exact_membership_registration_energy(memberships, priors)
    actual_energy = vectorized_toads_registration_energy(
        memberships,
        atlas,
        image_resolutions_mm=image_res,
        template_resolutions_mm=template_res,
        transform=transform,
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        scale=scale,
        image_center_xyz=image_center,
        template_center_xyz=template_center,
    )
    assert actual_energy == reference_energy


def test_scaled_registration_energy_and_coefficients_match_java_harness() -> None:
    fixture_dir = Path(__file__).parent / "data" / "synthetic" / "toads_scaled_registration_kernels"
    membership_image = nib.load(fixture_dir / "memberships.nii.gz")
    atlas_image = nib.load(fixture_dir / "atlas_shapes.nii.gz")
    assert nib.aff2axcodes(membership_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(atlas_image.affine) == ("R", "A", "S")
    memberships = np.asanyarray(membership_image.dataobj, dtype=np.float32)
    atlas = np.asanyarray(atlas_image.dataobj, dtype=np.float32)
    x, y, z, k = np.indices(memberships.shape, dtype=np.int32)
    np.testing.assert_array_equal(
        memberships,
        (((7 * x + 5 * y + 3 * z + 2 * k) % 17) / np.float32(16)).astype(np.float32),
    )
    x, y, z, k = np.indices(atlas.shape, dtype=np.int32)
    atlas_reference = (((2 * x + 5 * y + 3 * z + 4 * k) % 11) / np.float32(10)).astype(np.float32)
    atlas_reference[x >= 8] = 0.0
    np.testing.assert_array_equal(
        atlas,
        atlas_reference,
    )
    harness = Path(__file__).parent / "java_reference" / "ToadsScaledRegistrationKernelsReferenceHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == (
        "16347a6c4f2e1719c3d35368372a32404392e2959f020d1a01959eda2caac3f8"
    )
    reference = json.loads((fixture_dir / "java_reference.json").read_text())
    levels = build_toads_registration_pyramid(memberships, atlas, levels=3)
    assert [level.scale for level in levels] == [4, 2]
    image_center = tuple(size / 2.0 for size in memberships.shape[:3])
    template_center = tuple(size / 2.0 for size in atlas.shape[:3])
    transform = (0.14, 0.07, -0.05, 0.2, -0.1, 0.3)
    image_res = (0.8, 1.2, 1.5)
    template_res = (0.9, 1.1, 1.3)
    max_absolute_difference = 0.0
    for level in levels:
        expected = reference[f"scale_{level.scale}"]
        dimensions = level.memberships.shape[:3]
        low = np.zeros((3, 3), dtype=np.int64)
        high = np.broadcast_to(np.asarray(dimensions, dtype=np.int64), low.shape).copy()
        actual_energy = vectorized_toads_registration_energy(
            level.memberships,
            level.atlas_shapes,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            scale=level.scale,
            image_center_xyz=image_center,
            template_center_xyz=template_center,
        )
        actual = vectorized_toads_registration_coefficients(
            level.memberships,
            level.atlas_shapes,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            scale=level.scale,
            image_center_xyz=image_center,
            template_center_xyz=template_center,
        )
        scale_energy_cost_norm = np.asarray(expected["scale_energy_cost_norm"], np.float32)
        assert scale_energy_cost_norm[0] == np.float32(level.scale)
        for actual_value, reference_value, name in (
            (actual_energy, scale_energy_cost_norm[1], "energy"),
            (actual.cost, scale_energy_cost_norm[2], "cost"),
            (actual.norm, scale_energy_cost_norm[3], "norm"),
        ):
            difference = float(abs(actual_value - reference_value))
            max_absolute_difference = max(max_absolute_difference, difference)
            np.testing.assert_allclose(
                actual_value, reference_value, rtol=2e-7, atol=3e-8,
                err_msg=f"scale {level.scale} {name}",
            )
        for name in ("gradient", "hessian", "norm_gradient", "norm_hessian"):
            expected_values = np.asarray(expected[name], np.float32)
            difference = float(np.max(np.abs(getattr(actual, name) - expected_values)))
            max_absolute_difference = max(max_absolute_difference, difference)
            np.testing.assert_allclose(
                getattr(actual, name), expected_values, rtol=1e-5, atol=3e-8,
                err_msg=f"scale {level.scale} {name}",
            )
        positive_zero = np.asarray(expected["positive_zero_counts"], dtype=np.int64)
        assert positive_zero[0] > 0
        assert positive_zero[1] > 0
        assert positive_zero.sum() == int(np.prod(dimensions)) * 2
        scalar_positive = 0
        scalar_zero = 0
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    for class_index in (1, 2):
                        sample = sample_toad_shape_prior(
                            level.atlas_shapes[..., class_index],
                            (x, y, z),
                            image_dimensions_xyz=dimensions,
                            image_resolutions_mm=image_res,
                            template_resolutions_mm=template_res,
                            transform=transform,
                            class_index=class_index,
                            bbox_low_xyz=(0, 0, 0),
                            bbox_high_xyz=dimensions,
                            shape_slope=1.0,
                            scale=level.scale,
                            image_center_xyz=image_center,
                            template_center_xyz=template_center,
                        )
                        if sample.value > 0:
                            scalar_positive += 1
                        else:
                            scalar_zero += 1
        assert (scalar_positive, scalar_zero) == tuple(int(value) for value in positive_zero)
    assert max_absolute_difference <= 1.5e-7, max_absolute_difference
