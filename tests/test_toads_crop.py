from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_crop import (
    crop_toads_images,
    find_toads_crop_bounds,
    toads_rotation_matrix,
    transform_toads_topology,
    uncrop_toads_images,
)


FIXTURE = Path(__file__).parent / "data/synthetic/toads_crop_geometry/java_reference.json"


def _phantom():
    case = json.loads(FIXTURE.read_text())
    nx, ny, nz = case["dims_image"]
    tx, ty, tz = case["dims_template"]
    images = np.zeros((1, nx, ny, nz), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if 1 <= x <= 5 and 2 <= y <= 6 and 1 <= z <= 4 and (x + 2 * y + 3 * z) % 5:
                    images[0, x, y, z] = np.float32((x * 17 + y * 7 + z * 3) % 13) / np.float32(10)
    topology = np.ones((tx, ty, tz), dtype=np.uint8)
    for x in range(1, 4):
        for y in range(2, 5):
            for z in range(1, 4):
                if x in (1, 3) or y in (2, 4) or z in (1, 3):
                    topology[x, y, z] = 2
    topology[3:7, 3:7, 2:6] = 3
    transform = np.asarray(case["initial_transform"], dtype=np.float32)
    return case, images, topology, transform


def test_toads_crop_and_topology_match_source_formula_reference() -> None:
    case, images, topology, transform = _phantom()
    bounds = find_toads_crop_bounds(
        images,
        topology,
        transform,
        image_resolution=tuple(case["image_resolution"]),
        template_resolution=tuple(case["template_resolution"]),
        orientation_signs=tuple(case["image_orientation_signs"]),
        robust_ranges=np.asarray([[0.0, 1.0]], dtype=np.float32),
        background_ratio=0.2,
        background_label=1,
    )
    assert list(bounds) == case["bounds_inclusive"]

    cropped = crop_toads_images(
        images,
        bounds,
        case["padding"],
        transform,
        image_resolution=tuple(case["image_resolution"]),
        orientation_signs=tuple(case["image_orientation_signs"]),
    )
    assert list(cropped.images.shape[1:]) == case["cropped_dims"]
    np.testing.assert_array_equal(
        cropped.images[0], np.asarray(case["cropped_image"], dtype=np.float32)
    )
    np.testing.assert_allclose(cropped.transform, case["cropped_transform"], rtol=0, atol=2e-7)

    transformed = transform_toads_topology(
        topology,
        bounds,
        case["padding"],
        cropped.transform,
        crop_shape=tuple(case["cropped_dims"]),
        image_resolution=tuple(case["image_resolution"]),
        template_resolution=tuple(case["template_resolution"]),
        orientation_signs=tuple(case["image_orientation_signs"]),
        background_label=1,
    )
    np.testing.assert_array_equal(transformed, np.asarray(case["transformed_topology"]))

    uncropped, restored_transform = uncrop_toads_images(
        cropped.images,
        tuple(case["dims_image"]),
        bounds,
        case["padding"],
        cropped.transform,
        image_resolution=tuple(case["image_resolution"]),
        orientation_signs=tuple(case["image_orientation_signs"]),
    )
    np.testing.assert_array_equal(
        uncropped[0], np.asarray(case["uncropped_image"], dtype=np.float32)
    )
    np.testing.assert_allclose(restored_transform, case["uncropped_transform"], rtol=0, atol=2e-7)


@pytest.mark.parametrize("parameters", [np.asarray([np.nan, 0, 0]), np.zeros(4)])
def test_toads_rotation_rejects_malformed_parameters(parameters: np.ndarray) -> None:
    with pytest.raises(ValueError, match="three finite"):
        toads_rotation_matrix(parameters)
