from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_registration import (
    build_atlas_em_rigid_registration_pyramid,
    register_and_recrop_atlas_em,
    register_atlas_em_rigid_pyramid,
)
from pymacruise.atlas_em_preprocessing import (
    AtlasEMNormalizedImages,
    crop_atlas_em_images,
)
from pymacruise.toads_registration_sampling import register_toads_atlas_rigid_pyramid


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMRigidPyramidReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _field(salt: int) -> np.ndarray:
    values = np.empty((8, 8, 8, 3), dtype=np.float32)
    for x in range(8):
        for y in range(8):
            for z in range(8):
                for klass in range(3):
                    values[x, y, z, klass] = np.float32(
                        np.float32((x * 17 + y * 11 + z * 7 + klass * 13 + salt) % 97)
                        / np.float32(32.0)
                    )
    return values


def _bits(values: np.ndarray) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).ravel(order="C").view(np.uint32)]


def test_atlas_em_rigid_registration_pyramid_matches_java_subsampling_excerpt(
    tmp_path: Path,
) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-sourcepath", str(UPSTREAM), "-d", str(classes), str(HARNESS)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [java, "-cp", str(classes), "AtlasEMRigidPyramidReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    reference = json.loads(result.stdout)
    assert reference["demon_atlas_source_sha256"] == (
        "96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a"
    )
    assert reference["image_functions_sha256"] == (
        "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892"
    )

    actual = build_atlas_em_rigid_registration_pyramid(
        _field(3), _field(19), levels=3
    )
    assert [level.scale for level in actual] == [level["scale"] for level in reference["levels"]]
    for level, expected in zip(actual, reference["levels"], strict=True):
        assert list(level.memberships.shape[:3]) == expected["dims"]
        assert _bits(level.memberships) == expected["membership_bits"]
        assert _bits(level.atlas_shapes) == expected["shape_bits"]


def test_atlas_em_rigid_registration_selects_demon_source_profile() -> None:
    memberships = _field(3)
    memberships /= memberships.sum(axis=-1, keepdims=True, dtype=np.float32)
    shapes = _field(19)
    kwargs = {
        "levels": 2,
        "image_resolutions_mm": (1.0, 1.2, 1.5),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "initial_transform": (0.01, -0.02, 0.015, 0.0, 0.0, 0.0),
    }

    actual = register_atlas_em_rigid_pyramid(
        memberships, shapes, registration_iterations=2, **kwargs
    )
    reference = register_toads_atlas_rigid_pyramid(
        memberships,
        shapes,
        **kwargs,
        max_iterations=2,
        subsample_factor=2,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
        source_faithful=False,
        allow_nonfinite_trials=True,
    )
    assert [level.scale for level in actual.levels] == [2]
    assert np.array_equal(actual.transform, reference.transform)
    assert len(actual.levels[0].result.iterations) == len(reference.levels[0].result.iterations)
    for observed, expected in zip(
        actual.levels[0].result.iterations,
        reference.levels[0].result.iterations,
        strict=True,
    ):
        assert np.array_equal(observed.transform, expected.transform)
        assert observed.energy == expected.energy
        assert observed.lambda_value == expected.lambda_value
        assert observed.relative_change == expected.relative_change


def test_atlas_em_registration_hands_transform_into_native_recrop() -> None:
    axis = np.arange(12, dtype=np.float32)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    image = np.zeros((12, 12, 12, 1), dtype=np.float32)
    foreground = (x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2 < 17
    image[..., 0][foreground] = 0.5 + 0.02 * x[foreground]
    normalized = AtlasEMNormalizedImages(image, np.array([0.0], np.float32), np.array([1.0], np.float32))
    crop = crop_atlas_em_images(
        image,
        (2, 2, 2),
        (9, 9, 9),
        border=2,
        signed_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    cx, cy, cz = np.meshgrid(np.arange(12), np.arange(12), np.arange(12), indexing="ij")
    raw = np.stack(
        [
            np.ones((12, 12, 12), np.float32),
            np.exp(-((cx - 4) ** 2 + (cy - 5) ** 2 + (cz - 6) ** 2) / 18).astype(np.float32),
            np.exp(-((cx - 7) ** 2 + (cy - 6) ** 2 + (cz - 5) ** 2) / 18).astype(np.float32),
        ],
        axis=-1,
    )
    memberships = raw / raw.sum(axis=-1, keepdims=True)
    shapes = np.array(raw, dtype=np.float32)

    result = register_and_recrop_atlas_em(
        crop,
        memberships,
        shapes,
        normalized,
        background_ratio=0.1,
        levels=2,
        template_resolutions_mm=(1.0, 1.0, 1.0),
        registration_iterations=1,
    )
    assert [level.scale for level in result.registration.levels] == [2]
    assert len(result.registration.levels[0].result.iterations) == 1
    assert result.recropped.original_dimensions_xyz == (12, 12, 12)
    assert result.recropped.images.shape[3] == 1
    assert np.isfinite(result.registration.transform).all()
    assert np.isfinite(result.recropped.images).all()
