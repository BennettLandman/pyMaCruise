from __future__ import annotations

import numpy as np

from pymacruise.toads_shapes import compute_toad_shape_bounds, sample_toad_shape_prior
from pymacruise.toads_shape_volume import transform_toads_shape_priors


def test_vectorized_shape_volume_matches_selected_pointwise_sampler() -> None:
    template_dims = (8, 7, 6)
    image_dims = (6, 5, 4)
    x, y, z = np.indices(template_dims, dtype=np.float32)
    shape_a = np.clip(0.17 + 0.031 * x + 0.027 * y + 0.019 * z, 0.0, 1.0)
    shape_b = np.clip(0.83 - 0.021 * x + 0.013 * y - 0.024 * z, 0.0, 1.0)
    shapes = np.stack((shape_a, shape_b), axis=-1).astype(np.float32)
    image_spacing = (-0.8, -1.2, 1.5)
    template_spacing = (0.9, 1.1, 1.3)
    transform = (0.13, -0.07, 0.04, 0.35, -0.2, 0.15)
    bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., klass],
            image_dimensions_xyz=image_dims,
            image_resolutions_mm=image_spacing,
            template_resolutions_mm=template_spacing,
            transform=transform,
        )
        for klass in range(shapes.shape[-1])
    )
    actual = transform_toads_shape_priors(
        shapes,
        image_dimensions_xyz=image_dims,
        image_resolutions_mm=image_spacing,
        template_resolutions_mm=template_spacing,
        transform=transform,
        class_bounds=bounds,
    )

    expected = np.empty_like(actual)
    for x_index in range(image_dims[0]):
        for y_index in range(image_dims[1]):
            for z_index in range(image_dims[2]):
                for klass in range(shapes.shape[-1]):
                    expected[x_index, y_index, z_index, klass] = sample_toad_shape_prior(
                        shapes[..., klass],
                        (x_index, y_index, z_index),
                        image_dimensions_xyz=image_dims,
                        image_resolutions_mm=image_spacing,
                        template_resolutions_mm=template_spacing,
                        transform=transform,
                        class_index=klass,
                        bbox_low_xyz=bounds[klass].lower_xyz,
                        bbox_high_xyz=bounds[klass].upper_xyz,
                    ).value
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=3e-7)


def test_shape_volume_preserves_source_default_values_for_empty_support() -> None:
    shapes = np.zeros((4, 4, 4, 2), dtype=np.float32)
    dims = (3, 3, 3)
    bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., klass],
            image_dimensions_xyz=dims,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        )
        for klass in range(2)
    )
    result = transform_toads_shape_priors(
        shapes,
        image_dimensions_xyz=dims,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        class_bounds=bounds,
    )
    np.testing.assert_array_equal(result[..., 0], np.ones(dims, dtype=np.float32))
    np.testing.assert_array_equal(result[..., 1], np.zeros(dims, dtype=np.float32))
