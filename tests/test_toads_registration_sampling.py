from __future__ import annotations

import numpy as np
import pytest

from pymacruise.toads_registration import (
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_registration_sampling import (
    rigid_registration_coordinate_jacobians,
    sample_toads_registration_fields,
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
