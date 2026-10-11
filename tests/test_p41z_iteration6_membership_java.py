from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from compare_p41z_iteration6_membership_java import (
    HARNESS,
    _parse_java,
    _save_binary_input,
    compare,
)
from pymacruise.lesion_next_iteration_membership import run_lesion_toads_next_gm_membership_update
from pymacruise.lesion_next_iteration_statistics import (
    LesionNextIterationStatistics,
    LesionToadsNextStatisticsRefresh,
)
from pymacruise.lesion_statistics import LesionToadsAuxiliaryCentroids, LesionToadsVariances


def test_java_membership_energy_carries_from_active_lesion_to_skipped_support() -> None:
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for selected-source comparison")
    spatial = (4, 4, 3)
    images = np.empty((*spatial, 2), dtype=np.float32)
    images[..., 0] = np.float32(0.6)
    images[..., 1] = np.float32(0.8)
    k = np.arange(13, dtype=np.float32)
    centroids = np.stack((np.float32(0.2) + np.float32(0.03) * k,
                          np.float32(0.3) + np.float32(0.02) * k))
    stddev = np.full((2, 13), np.float32(0.25), dtype=np.float32)
    stats = LesionNextIterationStatistics(
        LesionToadsNextStatisticsRefresh(
            LesionToadsAuxiliaryCentroids(
                centroids, np.asarray([0.45, 0.7], dtype=np.float32), None,
            ),
            LesionToadsVariances(
                stddev, np.asarray([0.2, 0.3], dtype=np.float32),
                np.zeros(2, dtype=np.float32),
            ),
            np.asarray([0.5, 0.5], dtype=np.float32),
        ),
        ("synthetic-statistics",),
    )
    shape = np.ones((*spatial, 13), dtype=np.float32)
    shape[..., 10] = 0.0
    shape[1, 1, 1, 10] = np.float32(0.8)
    memberships = np.full((*spatial, 13), np.float32(0.02), dtype=np.float32)
    memberships[2, 1, 1, 10] = np.float32(0.15)
    lesion = np.full(spatial, np.float32(0.1), dtype=np.float32)
    zero = np.zeros(spatial, dtype=np.float32)
    arrays = {
        "images": images,
        "shape_values": shape,
        "direct_relations": np.ones((*spatial, 13), dtype=np.bool_),
        "class_memberships": memberships,
        "lesion": lesion,
        "black": zero.copy(),
        "outlier": zero.copy(),
        "lesion_relation": zero.copy(),
        "distance_maps": np.full((*spatial, 3), -1, dtype=np.int16),
    }
    with tempfile.TemporaryDirectory(prefix="p41z-lesion-energy-carry-") as temporary:
        work = Path(temporary)
        input_path, output_path = work / "input.bin", work / "output.bin"
        _save_binary_input(input_path, arrays, stats)
        classes = work / "classes"
        classes.mkdir()
        subprocess.run(["javac", "-d", str(classes), str(HARNESS)], check=True, capture_output=True)
        subprocess.run(
            ["java", "-Xmx512m", "-cp", str(classes), "LesionNextIterationMembershipReferenceHarness",
             str(input_path), str(output_path)],
            check=True, capture_output=True,
        )
        expected = _parse_java(output_path, spatial)

    actual = run_lesion_toads_next_gm_membership_update(
        stats, images, shape, np.ones(13, dtype=np.float32), np.ones(13, dtype=np.float32),
        arrays["direct_relations"], memberships, lesion, zero, zero, zero,
        arrays["distance_maps"][..., 0], arrays["distance_maps"][..., 1], arrays["distance_maps"][..., 2],
        previous_energy=0.0, outlier_energy=0.1, lesion_power_exponent=1.2,
    ).membership
    first, skipped = (1, 1, 1), (2, 1, 1)
    assert expected["lesion_energy"][first].view(np.uint32) != 0
    assert expected["lesion_energy"][skipped].view(np.uint32) == expected["lesion_energy"][first].view(np.uint32)
    assert expected["class_memberships"][skipped + (10,)] == np.float32(0.0)
    assert actual.class_memberships[skipped + (10,)] == np.float32(0.0)
    assert actual.lesion_energy[skipped].view(np.uint32) == actual.lesion_energy[first].view(np.uint32)
    assert actual.class_energies[skipped + (10,)].view(np.uint32) == actual.class_energies[first + (10,)].view(np.uint32)
    # The Java harness writes its accumulated objective after this voxel. The
    # compact fixture checks that it is finite and nonzero alongside the
    # carried per-voxel energy; full bitwise objective parity is exercised by
    # the opt-in fixed P41z case below.
    assert expected["lesion"][skipped] > np.float32(0.0)
    assert np.isfinite(expected["running_objective"][skipped])
    assert expected["running_objective"][skipped] > np.float32(0.0)


@pytest.mark.skipif(
    os.environ.get("PYMACRUISE_RUN_P41Z_ITER6_MEMBERSHIP") != "1",
    reason="set PYMACRUISE_RUN_P41Z_ITER6_MEMBERSHIP=1 for full-volume Java/Python comparison",
)
def test_p41z_iteration6_held_fixed_membership_matches_selected_java() -> None:
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for selected-source comparison")
    report = compare()
    assert report["spatial_shape_xyz"] == [71, 89, 67]
    assert report["statistics_exact_first"] is True
    assert report["nan_counts"]["java"] == report["nan_counts"]["python"]
    assert report["all_exact"], report
