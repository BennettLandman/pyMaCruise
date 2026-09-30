from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_shapes import sample_toad_shape_prior
from pymacruise.toads_registration import rigid_image_to_template_jacobian


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_shape_sampling"
    / "java_reference.json"
)
JAVA_SOURCE_HASHES = {
    "methods/ParametricTransform.java": "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac",
    "structures/RotationMatrix.java": "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0",
    "libraries/ImageFunctions.java": "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892",
    "utilities/Numerics.java": "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb",
}
SCALED_RIGID_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_scaled_rigid"
    / "java_reference.json"
)
SCALED_RIGID_HARNESS_SHA256 = "34a30727b2c61ef7b8e7bcaa8cd633605b6e989e7deeda91827f8e9877fc8eca"


def test_scaled_rigid_coordinates_and_jacobians_match_java_harness() -> None:
    fixture = json.loads(SCALED_RIGID_FIXTURE.read_text())
    harness = Path(__file__).parent / "java_reference" / "ToadsScaledRigidReferenceHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == SCALED_RIGID_HARNESS_SHA256
    shape = np.ones((11, 10, 9), dtype=np.float32)
    for case in fixture["cases"]:
        scale = case["scale"]
        point = tuple(fixture["point"])
        sample = sample_toad_shape_prior(
            shape,
            point,
            image_dimensions_xyz=tuple(fixture["image_dimensions"]),
            image_resolutions_mm=tuple(fixture["image_resolutions"]),
            template_resolutions_mm=tuple(fixture["template_resolutions"]),
            transform=tuple(fixture["transform"]),
            class_index=1,
            bbox_low_xyz=(0, 0, 0),
            bbox_high_xyz=tuple(size // int(scale) for size in fixture["image_dimensions"]),
            shape_slope=1.0,
            scale=scale,
            image_center_xyz=tuple(fixture["image_center"]),
            template_center_xyz=tuple(fixture["template_center"]),
        )
        np.testing.assert_array_equal(
            np.asarray(sample.template_xyz, dtype=np.float32),
            np.asarray(case["mapped"], dtype=np.float32),
        )
        jacobian = rigid_image_to_template_jacobian(
            point,
            fixture["transform"],
            fixture["image_center"],
            fixture["image_resolutions"],
            fixture["template_resolutions"],
            scale=scale,
        )
        np.testing.assert_array_equal(jacobian, np.asarray(case["jacobian"], np.float32))


@pytest.mark.parametrize(
    ("scale", "image_center", "template_center", "message"),
    [
        (0.0, (4.5, 4.0, 3.5), (5.5, 5.0, 4.5), "scale must be positive"),
        (2.0, (4.5, 4.0), (5.5, 5.0, 4.5), "image_center_xyz must contain three"),
        (2.0, (4.5, 4.0, 3.5), (5.5, 5.0, float("nan")), "template_center_xyz must contain three"),
    ],
)
def test_scaled_shape_sampler_rejects_invalid_scale_and_centers(
    scale: float,
    image_center: tuple[float, ...],
    template_center: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        sample_toad_shape_prior(
            np.ones((11, 10, 9), dtype=np.float32),
            (1, 1, 0),
            image_dimensions_xyz=(9, 8, 7),
            image_resolutions_mm=(0.8, 1.2, 2.0),
            template_resolutions_mm=(0.9, 1.1, 1.3),
            transform=(0.14, 0.07, -0.05, 1.0, 2.0, 3.0),
            class_index=1,
            bbox_low_xyz=(0, 0, 0),
            bbox_high_xyz=(4, 4, 4),
            scale=scale,
            image_center_xyz=image_center,
            template_center_xyz=template_center,
        )


def _shape_prior(template_dims: tuple[int, int, int]) -> np.ndarray:
    x, y, z = np.indices(template_dims, dtype=np.float32)
    return (0.2 + 0.02 * x + 0.025 * y + 0.03 * z + 0.005 * x * y).astype(np.float32)


def test_toads_shape_sampler_matches_selected_java_transform_and_interpolation() -> None:
    reference = json.loads(FIXTURE.read_text())
    shape = _shape_prior(tuple(reference["template_dims"]))
    for case in reference["cases"]:
        actual = sample_toad_shape_prior(
            shape,
            tuple(case["point"]),
            image_dimensions_xyz=tuple(reference["image_dims"]),
            image_resolutions_mm=tuple(reference["image_resolutions"]),
            template_resolutions_mm=tuple(reference["template_resolutions"]),
            transform=tuple(case["transform"]),
            class_index=case["class"],
            bbox_low_xyz=tuple(case["bbox_low"]),
            bbox_high_xyz=tuple(case["bbox_high"]),
            shape_slope=case["slope"],
        )
        assert actual.value == np.float32(case["value"]), case["name"]
        if case["mapped"] is None:
            assert actual.template_xyz is None
        else:
            np.testing.assert_array_equal(
                np.asarray(actual.template_xyz, dtype=np.float32),
                np.asarray(case["mapped"], dtype=np.float32),
                err_msg=case["name"],
            )


@pytest.mark.parametrize(("relative_path", "sha256"), JAVA_SOURCE_HASHES.items())
def test_toads_shape_reference_sources_match_recorded_hashes(
    relative_path: str, sha256: str
) -> None:
    source = (
        Path(__file__).parent
        / "java_reference"
        / "upstream"
        / "edu"
        / "jhmi"
        / "rad"
        / "medic"
        / relative_path
    )
    assert hashlib.sha256(source.read_bytes()).hexdigest() == sha256


@pytest.mark.parametrize(
    ("shape", "point", "image_dims", "image_res", "template_res", "transform", "low", "high", "message"),
    [
        (np.ones((4, 4)), (1, 1, 1), (4, 4, 4), (1, 1, 1), (1, 1, 1), (0, 0, 0, 0, 0, 0), (0, 0, 0), (4, 4, 4), "3D"),
        (np.ones((4, 4, 4)), (4, 1, 1), (4, 4, 4), (1, 1, 1), (1, 1, 1), (0, 0, 0, 0, 0, 0), (0, 0, 0), (4, 4, 4), "inside"),
        (np.ones((4, 4, 4)), (1, 1, 1), (4, 4, 4), (0, 1, 1), (1, 1, 1), (0, 0, 0, 0, 0, 0), (0, 0, 0), (4, 4, 4), "nonzero signed"),
        (np.ones((4, 4, 4)), (1, 1, 1), (4, 4, 4), (1, 1, 1), (1, 1, 1), (0, 0, 0), (0, 0, 0), (4, 4, 4), "six finite"),
    ],
)
def test_toads_shape_sampler_rejects_invalid_inputs(
    shape: np.ndarray,
    point: tuple[int, int, int],
    image_dims: tuple[int, int, int],
    image_res: tuple[float, float, float],
    template_res: tuple[float, float, float],
    transform: tuple[float, ...],
    low: tuple[int, int, int],
    high: tuple[int, int, int],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        sample_toad_shape_prior(
            shape,
            point,
            image_dimensions_xyz=image_dims,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=transform,
            class_index=1,
            bbox_low_xyz=low,
            bbox_high_xyz=high,
        )
