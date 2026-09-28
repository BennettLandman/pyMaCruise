from __future__ import annotations

import numpy as np
import pytest

from pymacruise.toads_registration import (
    exact_membership_registration_coefficients,
    exact_membership_registration_energy,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_registration_sampling import (
    register_toads_atlas_rigid,
    rigid_registration_coordinate_jacobians,
    sample_toads_registration_fields,
    stream_toads_registration_coefficients,
    stream_toads_registration_energy,
)
from pymacruise.toads_shapes import sample_toad_shape_prior


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
