from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_registration import register_atlas_em_rigid_single_scale


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
    assert 1 <= len(result.iterations) <= iterations
