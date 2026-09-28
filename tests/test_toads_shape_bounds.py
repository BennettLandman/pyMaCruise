from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

from pymacruise.toads_shapes import compute_toad_shape_bounds, sample_toad_shape_prior


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_shape_bounds"
    / "java_reference.json"
)
ATLAS_SOURCE_SHA256 = "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2"
TRANSFORM_SOURCE_SHA256 = "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac"
FOLDED_FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_shape_bounds_folded"
    / "java_reference.json"
)


def _prior(kind: str, dims: tuple[int, int, int]) -> np.ndarray:
    result = np.empty(dims, dtype=np.float32)
    for x in range(dims[0]):
        for y in range(dims[1]):
            for z in range(dims[2]):
                if kind == "zero":
                    value = np.float32(0.0)
                elif kind == "positive":
                    value = np.float32(0.25)
                elif 1 <= x <= 4 and 1 <= y <= 3 and 1 <= z <= 2:
                    value = np.float32(
                        np.float32(
                            np.float32(np.float32(0.1) + np.float32(0.03) * x)
                            + np.float32(0.02) * y
                        )
                        + np.float32(0.01) * z
                    )
                else:
                    value = np.float32(-0.2)
                result[x, y, z] = value
    return result


def test_transformed_shape_support_bounds_match_java_scan_excerpt() -> None:
    reference = json.loads(FIXTURE.read_text())
    assert reference["source_sha256"] == ATLAS_SOURCE_SHA256
    assert reference["parametric_transform_sha256"] == TRANSFORM_SOURCE_SHA256
    image_dims = tuple(reference["image_dims"])
    template_dims = tuple(reference["template_dims"])
    for case in reference["cases"]:
        actual = compute_toad_shape_bounds(
            _prior(case["shape_kind"], template_dims),
            image_dimensions_xyz=image_dims,
            image_resolutions_mm=tuple(reference["image_resolutions"]),
            template_resolutions_mm=tuple(reference["template_resolutions"]),
            transform=tuple(case["transform"]),
            record_support=True,
        )
        assert list(actual.lower_xyz) == case["lower"], case["name"]
        assert list(actual.upper_xyz) == case["upper"], case["name"]
        assert actual.positive_support_xyz is not None
        assert [list(point) for point in actual.positive_support_xyz] == case["support"], (
            case["name"]
        )
        assert actual.support_count == case["support_count"]


def test_folded_prior_support_bounds_refresh_match_java_scan() -> None:
    reference = json.loads(FOLDED_FIXTURE.read_text())
    assert reference["source_sha256"] == ATLAS_SOURCE_SHA256
    assert reference["parametric_transform_sha256"] == TRANSFORM_SOURCE_SHA256
    prior_path = FOLDED_FIXTURE.parent / "folded_shape_prior.nii.gz"
    image = nib.load(str(prior_path))
    assert image.shape == tuple(reference["template_dims"])
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    prior = np.asarray(image.dataobj, dtype=np.float32)
    assert np.isfinite(prior).all()
    assert np.min(prior) == np.float32(0.0)
    assert np.max(prior) > np.float32(0.99)

    observed_bounds = []
    for case in reference["cases"]:
        actual = compute_toad_shape_bounds(
            prior,
            image_dimensions_xyz=tuple(reference["image_dims"]),
            image_resolutions_mm=tuple(reference["resolutions"]),
            template_resolutions_mm=tuple(reference["resolutions"]),
            transform=tuple(case["transform"]),
        )
        assert list(actual.lower_xyz) == case["lower"], case["name"]
        assert list(actual.upper_xyz) == case["upper"], case["name"]
        assert actual.support_count == case["support_count"], case["name"]
        observed_bounds.append((actual.lower_xyz, actual.upper_xyz))

    assert len(set(observed_bounds)) == len(observed_bounds)


def test_shape_sampler_preserves_source_maximum_as_exclusive_bound() -> None:
    reference = json.loads(FIXTURE.read_text())
    case = next(c for c in reference["cases"] if c["name"] == "all_positive_touches_all_faces")
    bounds = compute_toad_shape_bounds(
        _prior(case["shape_kind"], tuple(reference["template_dims"])),
        image_dimensions_xyz=tuple(reference["image_dims"]),
        image_resolutions_mm=tuple(reference["image_resolutions"]),
        template_resolutions_mm=tuple(reference["template_resolutions"]),
        transform=tuple(case["transform"]),
    )
    assert bounds.positive_support_xyz is None
    assert bounds.support_count == case["support_count"]
    kwargs = {
        "image_dimensions_xyz": tuple(reference["image_dims"]),
        "image_resolutions_mm": tuple(reference["image_resolutions"]),
        "template_resolutions_mm": tuple(reference["template_resolutions"]),
        "transform": tuple(case["transform"]),
        "class_index": 1,
        "bbox_low_xyz": bounds.lower_xyz,
        "bbox_high_xyz": bounds.upper_xyz,
    }
    included = sample_toad_shape_prior(
        _prior(case["shape_kind"], tuple(reference["template_dims"])),
        bounds.lower_xyz,
        **kwargs,
    )
    excluded_max_plane = sample_toad_shape_prior(
        _prior(case["shape_kind"], tuple(reference["template_dims"])),
        bounds.upper_xyz,
        **kwargs,
    )
    assert included.template_xyz is not None
    assert excluded_max_plane.template_xyz is None
    assert excluded_max_plane.value == np.float32(0.0)


def test_shape_sampler_accepts_the_source_empty_support_sentinel() -> None:
    reference = json.loads(FIXTURE.read_text())
    case = next(c for c in reference["cases"] if c["name"] == "no_positive_support")
    bounds = compute_toad_shape_bounds(
        _prior(case["shape_kind"], tuple(reference["template_dims"])),
        image_dimensions_xyz=tuple(reference["image_dims"]),
        image_resolutions_mm=tuple(reference["image_resolutions"]),
        template_resolutions_mm=tuple(reference["template_resolutions"]),
        transform=tuple(case["transform"]),
    )
    result = sample_toad_shape_prior(
        _prior(case["shape_kind"], tuple(reference["template_dims"])),
        (0, 0, 0),
        image_dimensions_xyz=tuple(reference["image_dims"]),
        image_resolutions_mm=tuple(reference["image_resolutions"]),
        template_resolutions_mm=tuple(reference["template_resolutions"]),
        transform=tuple(case["transform"]),
        class_index=0,
        bbox_low_xyz=bounds.lower_xyz,
        bbox_high_xyz=bounds.upper_xyz,
    )
    assert result.value == np.float32(1.0)
    assert result.template_xyz is None


@pytest.mark.parametrize(
    ("shape", "dims", "message"),
    [
        (np.zeros((3, 3), dtype=np.float32), (3, 3, 3), "3D"),
        (np.zeros((3, 3, 3), dtype=np.float32), (3, 3), "three integers"),
    ],
)
def test_shape_bounds_reject_invalid_arrays_or_grid(
    shape: np.ndarray, dims: tuple[int, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        compute_toad_shape_bounds(
            shape,
            image_dimensions_xyz=dims,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        )
