from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_registration import (
    register_atlas_em_rigid_pyramid,
    register_atlas_em_rigid_single_scale,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_RUNNER = PROJECT_ROOT / "tests/java_reference/run_atlas_em_demon_rigid_registration_reference.sh"


@pytest.mark.parametrize("iterations", (1, 3))
def test_atlas_em_rigid_update_matches_selected_demon_class(
    tmp_path: Path, iterations: int
) -> None:
    if shutil.which("javac") is None or shutil.which("java") is None:
        pytest.skip("JDK is needed to execute the selected AtlasEM Java class")

    dimensions = (13, 11, 9)
    x, y, z = np.indices(dimensions, dtype=np.float32)
    atlas = np.zeros((*dimensions, 2), dtype=np.float32)
    atlas[..., 0] = np.float32(0.25)
    atlas[..., 1] = np.exp(
        -(((x - 6.0) / 3.2) ** 2 + ((y - 5.0) / 2.5) ** 2 + ((z - 4.0) / 2.1) ** 2)
    ).astype(np.float32)
    target = np.exp(
        -(((x - 6.7) / 3.2) ** 2 + ((y - 4.4) / 2.5) ** 2 + ((z - 4.5) / 2.1) ** 2)
    ).astype(np.float32)
    memberships = np.stack((np.float32(1.0) - target, target), axis=-1).astype(np.float32)

    # Java consumes class-major shape volumes and voxel-major membership classes.
    np.moveaxis(atlas, -1, 0).astype("<f4").tofile(tmp_path / "shapes.f32le")
    memberships.astype("<f4").tofile(tmp_path / "memberships.f32le")
    output = subprocess.check_output(
        [
            str(REFERENCE_RUNNER),
            *(str(value) for value in dimensions),
            str(iterations),
            "1",
            "single",
            str(tmp_path),
        ],
        cwd=PROJECT_ROOT,
        text=True,
        stderr=subprocess.DEVNULL,
    ).strip()
    java_transform = np.asarray(
        [int(value) for value in output.split(",")], dtype=np.int32
    ).view(np.float32)

    result = register_atlas_em_rigid_single_scale(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.2, 1.5),
        template_resolutions_mm=(0.8, 1.0, 1.25),
        initial_transform=(0.0,) * 6,
        registration_iterations=iterations,
        registered_classes=np.array([False, True]),
    )
    assert result.transform.dtype == np.float32
    np.testing.assert_array_equal(result.transform, java_transform)
    assert np.any(result.transform != np.zeros(6, dtype=np.float32))
    assert 1 <= len(result.iterations) <= iterations


@pytest.mark.parametrize("iterations", (1, 2, 5))
def test_atlas_em_pyramid_matches_selected_demon_class(
    tmp_path: Path, iterations: int
) -> None:
    if shutil.which("javac") is None or shutil.which("java") is None:
        pytest.skip("JDK is needed to execute the selected AtlasEM Java class")

    dimensions = (16, 16, 16)
    x, y, z = np.indices(dimensions, dtype=np.float32)
    atlas = np.zeros((*dimensions, 2), dtype=np.float32)
    atlas[..., 0] = np.float32(0.25)
    atlas[..., 1] = np.exp(
        -(((x - 7.5) / 4.0) ** 2 + ((y - 7.0) / 3.5) ** 2 + ((z - 8.0) / 3.0) ** 2)
    ).astype(np.float32)
    target = np.exp(
        -(((x - 8.0) / 4.0) ** 2 + ((y - 6.5) / 3.5) ** 2 + ((z - 8.5) / 3.0) ** 2)
    ).astype(np.float32)
    memberships = np.stack((np.float32(1.0) - target, target), axis=-1).astype(np.float32)

    np.moveaxis(atlas, -1, 0).astype("<f4").tofile(tmp_path / "shapes.f32le")
    memberships.astype("<f4").tofile(tmp_path / "memberships.f32le")
    np.array([False, True], dtype=np.uint8).tofile(tmp_path / "registered.u8")
    output = subprocess.check_output(
        [
            str(REFERENCE_RUNNER),
            *(str(value) for value in dimensions),
            str(iterations),
            "3",
            "pyramid",
            str(tmp_path),
        ],
        cwd=PROJECT_ROOT,
        text=True,
        stderr=subprocess.DEVNULL,
    ).strip()
    java_transform = np.asarray(
        [int(value) for value in output.split(",")], dtype=np.int32
    ).view(np.float32)

    result = register_atlas_em_rigid_pyramid(
        memberships,
        atlas,
        levels=3,
        image_resolutions_mm=(1.0, 1.2, 1.5),
        template_resolutions_mm=(0.8, 1.0, 1.25),
        initial_transform=(0.0,) * 6,
        registration_iterations=iterations,
        registered_classes=np.array([False, True]),
    )
    np.testing.assert_array_equal(result.transform, java_transform)
    assert np.any(result.transform != np.zeros(6, dtype=np.float32))
    assert len(result.levels) == 2


@pytest.mark.parametrize(
    "registered_flags",
    ((False, True, False), (False, True, True)),
)
def test_atlas_em_pyramid_matches_selected_demon_for_selected_classes(
    tmp_path: Path,
    registered_flags: tuple[bool, bool, bool],
) -> None:
    if shutil.which("javac") is None or shutil.which("java") is None:
        pytest.skip("JDK is needed to execute the selected AtlasEM Java class")

    dimensions = (24, 24, 24)
    x, y, z = np.indices(dimensions, dtype=np.float32)
    shapes = np.zeros((*dimensions, 3), dtype=np.float32)
    shapes[..., 0] = np.float32(0.25)
    shapes[..., 1] = np.exp(-(((x - 8.0) / 3.2) ** 2 + ((y - 12.0) / 3.1) ** 2 + ((z - 11.0) / 3.4) ** 2))
    shapes[..., 2] = np.exp(-(((x - 15.0) / 3.1) ** 2 + ((y - 10.0) / 3.6) ** 2 + ((z - 14.0) / 2.9) ** 2))
    shapes = shapes.astype(np.float32)
    memberships = np.empty_like(shapes)
    memberships[..., 0] = np.float32(0.3)
    memberships[..., 1] = np.float32(0.4)
    memberships[..., 2] = np.float32(0.3)
    memberships[..., 1] += np.float32(0.1) * np.exp(-(((x - 9.0) / 5.0) ** 2 + ((y - 11.0) / 5.0) ** 2 + ((z - 12.0) / 5.0) ** 2))
    memberships[..., 2] -= memberships[..., 1] - np.float32(0.4)
    memberships = memberships.astype(np.float32)
    registered = np.asarray(registered_flags, dtype=np.uint8)

    np.moveaxis(shapes, -1, 0).astype("<f4").tofile(tmp_path / "shapes.f32le")
    memberships.astype("<f4").tofile(tmp_path / "memberships.f32le")
    registered.tofile(tmp_path / "registered.u8")
    output = subprocess.check_output(
        [
            str(REFERENCE_RUNNER),
            *(str(value) for value in dimensions),
            "3",
            "3",
            "pyramid",
            str(tmp_path),
        ],
        cwd=PROJECT_ROOT,
        text=True,
        stderr=subprocess.DEVNULL,
    ).strip()
    java_transform = np.asarray(
        [int(value) for value in output.split(",")], dtype=np.int32
    ).view(np.float32)
    result = register_atlas_em_rigid_pyramid(
        memberships,
        shapes,
        levels=3,
        image_resolutions_mm=(1.0, 1.2, 1.5),
        template_resolutions_mm=(0.8, 1.0, 1.25),
        initial_transform=(0.0,) * 6,
        registration_iterations=3,
        registered_classes=registered.astype(bool),
    )
    np.testing.assert_array_equal(result.transform, java_transform)
    assert np.any(result.transform != np.zeros(6, dtype=np.float32))
