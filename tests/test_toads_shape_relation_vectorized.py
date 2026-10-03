from __future__ import annotations

import numpy as np

from pymacruise.toads_math import toads_shape_relation_factors


def _scalar_reference(
    shapes: np.ndarray, labels: np.ndarray, prior_factor: float
) -> tuple[np.ndarray, np.ndarray]:
    """Original scalar control flow, retained here as an exact oracle."""
    factor = np.float32(prior_factor)
    out = np.empty(shapes.shape, dtype=np.float32)
    suppressed = np.zeros(shapes.shape, dtype=np.bool_)
    for x in range(shapes.shape[0]):
        for y in range(shapes.shape[1]):
            for z in range(shapes.shape[2]):
                for klass in range(shapes.shape[3]):
                    value = shapes[x, y, z, klass]
                    direct = bool(np.any(labels[x, y, z, :2] == klass))
                    neighbor = bool(np.any(labels[x, y, z, 2:] == klass))
                    if direct:
                        result = np.float32(1.0)
                    elif neighbor:
                        scaled = np.float32(np.float32(2.0) * factor)
                        scaled = np.float32(scaled * value)
                        result = np.float32(
                            max(np.float32(0.25), min(scaled, np.float32(1.0)))
                        )
                    elif value == np.float32(0.0):
                        result = np.float32(0.0)
                        suppressed[x, y, z, klass] = True
                    else:
                        scaled = np.float32(np.float32(2.0) * factor)
                        scaled = np.float32(scaled * value)
                        result = np.float32(min(np.float32(1.0), scaled))
                    out[x, y, z, klass] = result
    return out, suppressed


def test_vectorized_shape_relation_factors_are_bitwise_scalar_equivalent() -> None:
    rng = np.random.default_rng(348901)
    shapes = rng.uniform(0.0, 1.25, size=(13, 12, 11, 11)).astype(np.float32)
    # Exercise exact zero support, values around the clamps, and signed zero.
    shapes.flat[::17] = np.float32(0.0)
    shapes.flat[1::29] = np.nextafter(np.float32(0.25), np.float32(0.0))
    shapes.flat[2::31] = np.nextafter(np.float32(1.0), np.float32(2.0))
    shapes.flat[3::37] = np.float32(-0.0)
    labels = rng.integers(-1, 11, size=(13, 12, 11, 4), dtype=np.int16)
    # Force cases with direct/neighbor overlap and all-empty labels.
    labels[0, 0, 0] = (3, 3, 3, -1)
    labels[0, 0, 1] = (-1, -1, -1, -1)
    labels[0, 0, 2] = (2, 7, 4, 2)

    for prior_factor in (0.0, 0.125, 0.75, 1.0, -0.5):
        expected, expected_suppressed = _scalar_reference(shapes, labels, prior_factor)
        actual, actual_suppressed = toads_shape_relation_factors(
            shapes, labels, prior_factor
        )
        np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32))
        np.testing.assert_array_equal(actual_suppressed, expected_suppressed)


def test_vectorized_shape_relation_factors_preserve_empty_and_classless_edges() -> None:
    shapes = np.empty((0, 3, 2, 4), dtype=np.float32)
    labels = np.empty((0, 3, 2, 4), dtype=np.int16)
    actual, suppressed = toads_shape_relation_factors(shapes, labels, 0.5)
    assert actual.shape == shapes.shape and actual.dtype == np.float32
    assert suppressed.shape == shapes.shape and suppressed.dtype == np.bool_

    shapes = np.zeros((1, 1, 1, 4), dtype=np.float32)
    labels = np.full((1, 1, 1, 4), -1, dtype=np.int16)
    expected, expected_suppressed = _scalar_reference(shapes, labels, 0.5)
    actual, suppressed = toads_shape_relation_factors(shapes, labels, 0.5)
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(suppressed, expected_suppressed)

