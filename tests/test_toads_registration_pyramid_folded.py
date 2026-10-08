from __future__ import annotations

import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_registration_sampling import (
    register_toads_atlas_rigid_pyramid,
)


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_registration_pyramid_folded"
)


def test_folded_cortex_pyramid_matches_actual_vandycruise_source(
    tmp_path: Path,
) -> None:
    memberships_image = nib.load(FIXTURE / "memberships.nii.gz")
    shapes_image = nib.load(FIXTURE / "atlas-shape-priors.nii.gz")
    memberships = np.asarray(memberships_image.dataobj, dtype=np.float32)
    shapes = np.asarray(shapes_image.dataobj, dtype=np.float32)
    assert memberships.shape == shapes.shape == (48, 48, 48, 4)
    assert nib.aff2axcodes(memberships_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(shapes_image.affine) == ("R", "A", "S")

    memberships_path = tmp_path / "memberships-f32le.bin"
    shapes_path = tmp_path / "shapes-f32le.bin"
    memberships.tofile(memberships_path)
    shapes.tofile(shapes_path)
    dimensions = memberships.shape[:3]
    spacings = memberships_image.header.get_zooms()[:3]
    command = [
        str(
            Path(__file__).parent
            / "java_reference"
            / "run_toads_pyramid_actual_registration_reference.sh"
        ),
        str(memberships_path),
        str(shapes_path),
        *(str(value) for value in dimensions),
        *(str(value) for value in shapes.shape[:3]),
        str(memberships.shape[3]),
        *(str(float(value)) for value in spacings),
        *(str(float(value)) for value in spacings),
        "4",
        "8",
        *("0" for _ in range(6)),
    ]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    java_result = json.loads(
        next(line for line in completed.stdout.splitlines() if line.startswith("{"))
    )

    result = register_toads_atlas_rigid_pyramid(
        memberships,
        shapes,
        levels=4,
        image_resolutions_mm=tuple(float(value) for value in spacings),
        template_resolutions_mm=tuple(float(value) for value in spacings),
        initial_transform=(0.0,) * 6,
        registered_classes=np.asarray((False, True, True, True)),
        max_iterations=8,
    )
    np.testing.assert_array_equal(
        result.transform.view(np.int32),
        np.asarray(java_result["transform_bits"], dtype=np.int32),
    )

    source_levels = (4, 3, 2)
    for python_level, source_level in zip(result.levels, source_levels, strict=True):
        start = completed.stdout.index(f"level {source_level}\n")
        next_level = source_level - 1
        end_marker = f"level {next_level}\n"
        end = completed.stdout.find(end_marker, start + 1)
        section = completed.stdout[start:] if end < 0 else completed.stdout[start:end]
        assert section.count("trial transform:") == sum(
            item.rejected_trials + int(item.accepted)
            for item in python_level.result.iterations
        )
        assert section.count(" (E:") == sum(
            item.accepted for item in python_level.result.iterations
        )
        assert all(np.isfinite(item.transform).all() for item in python_level.result.iterations)
        accepted = [item for item in python_level.result.iterations if item.accepted]
        assert accepted

    finest_iterations = result.levels[-1].result.iterations
    finest_accepted = [item for item in finest_iterations if item.accepted][-1]
    assert finest_accepted.energy.view(np.int32) == np.float32(
        java_result["final_energy"]
    ).view(np.int32)
    # The Java object retains the lambda after the terminal rejected trial;
    # the accepted transform and its energy remain unchanged.
    assert finest_iterations[-1].lambda_value.view(np.int32) == np.int32(
        java_result["lambda_bits"]
    )
