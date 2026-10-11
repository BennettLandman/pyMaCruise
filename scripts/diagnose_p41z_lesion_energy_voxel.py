#!/usr/bin/env python3
"""Compare class-10 lesion-energy terms at fixed P41z voxel [1,28,29]."""
from __future__ import annotations

import json
import struct
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from compare_p41z_iteration6_membership_java import _stats_from_fixture
from pymacruise.lesion_membership import (
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
)
from pymacruise.lesion_similarity import compute_lesion_toads_similarity
from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies
from pymacruise.lesion_standard13_voxel import (
    STANDARD13_CLASS_SMOOTHING,
    STANDARD13_NAMES,
    STANDARD13_PRIOR_COEFFICIENT,
    STANDARD13_PRIOR_SCALE,
    STANDARD13_TYPES,
)
from pymacruise.toads_math import toads_powertable_lookup

ROOT = Path(__file__).resolve().parents[1]
JAVA = ROOT / "tests/java_reference/P41zSingleVoxelLesionEnergyDiagnostic.java"
TOADS = ROOT.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic"
POWER_TABLE = TOADS / "utilities/PowerTable.java"
COORD = (1, 28, 29)
FACES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def _f32bits(value: float | np.float32) -> str:
    return struct.pack(">f", np.float32(value)).hex()


def _write_input(path: Path, arrays: dict[str, np.ndarray], stats, case: dict) -> None:
    x, y, z = COORD
    memberships = arrays["class_memberships"]
    neighbors = np.asarray([memberships[x+dx, y+dy, z+dz] for dx, dy, dz in FACES], dtype=np.float32)
    lesion_faces = np.asarray([arrays["lesion"][x+dx, y+dy, z+dz] for dx, dy, dz in FACES], dtype=np.float32)
    black_faces = np.asarray([arrays["black"][x+dx, y+dy, z+dz] for dx, dy, dz in FACES], dtype=np.float32)
    settings = case["membership_settings"]
    values = [
        arrays["images"][x, y, z],
        stats.statistics.centroids.class_centroids,
        stats.statistics.variances.class_standard_deviation,
        stats.statistics.centroids.lesion_centroids,
        stats.statistics.variances.lesion_standard_deviation,
        stats.statistics.channel_weights,
        arrays["shape_values"][x, y, z],
        np.asarray(settings["shape_priors"], dtype=np.float32),
        np.asarray(settings["smoothing_priors"], dtype=np.float32),
        neighbors, lesion_faces, black_faces,
        np.asarray([settings["lesion_power_exponent"], float(STANDARD13_CLASS_SMOOTHING),
                    float(STANDARD13_PRIOR_COEFFICIENT)], dtype=np.float32),
    ]
    with path.open("wb") as stream:
        for value in values:
            stream.write(np.asarray(value, dtype=">f4").tobytes(order="C"))


def diagnose(java: str = "java", javac: str = "javac") -> dict:
    arrays, stats, case = _stats_from_fixture()
    x, y, z = COORD
    settings = case["membership_settings"]
    image = arrays["images"][x, y, z]
    shapes = arrays["shape_values"][x, y, z]
    class_neighbors = np.asarray([
        arrays["class_memberships"][x+dx, y+dy, z+dz] for dx, dy, dz in FACES
    ], dtype=np.float32)
    lesion_faces = np.asarray([arrays["lesion"][x+dx, y+dy, z+dz] for dx, dy, dz in FACES], dtype=np.float32)
    black_faces = np.asarray([arrays["black"][x+dx, y+dy, z+dz] for dx, dy, dz in FACES], dtype=np.float32)
    k = 10
    exponent = float(settings["lesion_power_exponent"])
    lesion_centroids = stats.statistics.centroids.lesion_centroids
    class_centroids = stats.statistics.centroids.class_centroids
    class_stddev = stats.statistics.variances.class_standard_deviation
    lesion_stddev = stats.statistics.variances.lesion_standard_deviation
    channel_weights = stats.statistics.channel_weights

    raw_result = compute_lesion_toads_unfielded_data_energies(
        image, class_centroids[:, k], lesion_centroids, np.zeros(2, dtype=np.float32),
        channel_weights, object_type="lesion", flair_channel=1,
        use_lesion_weight=bool(settings["use_lesion_weight"]),
        lesion_power_exponent=exponent, cached_lesion_power_exponent=exponent,
    )
    raw = np.float32(raw_result.lesion_energy)
    smooth_result = update_lesion_toads_smoothing_energies(
        class_neighbors, lesion_faces, black_faces, class_index=k, class_type="LESION",
        smoothing_prior=float(settings["smoothing_priors"][k]),
        smoothing=float(STANDARD13_CLASS_SMOOTHING), class_energy=0.0,
        lesion_energy=float(raw), black_hole_energy=0.0,
        class_membership_data_term=0.0, lesion_membership_data_term=float(raw),
        black_hole_membership_data_term=0.0, outmask_membership_data_term=0.0,
    )
    smooth_total = np.float32(smooth_result.lesion_energy)
    smooth_sum = np.float32(smooth_result.other_class_neighbor_squares
                            + smooth_result.current_class_neighbor_squares)
    smooth_sum = np.float32(smooth_sum + smooth_result.black_hole_neighbor_squares)
    smooth_increment = np.float32(
        np.float32(np.float32(np.float32(0.5) * np.float32(settings["smoothing_priors"][k]))
                   * np.float32(STANDARD13_CLASS_SMOOTHING)) * smooth_sum
    )
    similarities = compute_lesion_toads_similarity(
        class_centroids, class_stddev, lesion_centroids, lesion_stddev, channel_weights,
        STANDARD13_TYPES, STANDARD13_NAMES, prior_scale=float(STANDARD13_PRIOR_SCALE),
        t1_channel=0, black_hole_similarity=np.zeros(13, dtype=np.float32),
        finite_zero_variance=True,
    )
    prior_result = update_lesion_toads_prior_energies(
        np.asarray(settings["shape_priors"], dtype=np.float32),
        similarities.class_similarity[k], similarities.lesion_similarity,
        similarities.black_hole_similarity, shapes, class_index=k, class_type="LESION",
        prior_coefficient=float(STANDARD13_PRIOR_COEFFICIENT),
        class_membership_data_term=0.0,
        lesion_membership_data_term=float(smooth_result.lesion_membership_data_term),
        black_hole_membership_data_term=0.0, outmask_membership_data_term=0.0,
        class_energy=0.0, lesion_energy=float(smooth_total), black_hole_energy=0.0,
        has_lesion=True, lesion_class_index=10, lesion_class_type="LESION",
    )
    final_energy = np.float32(prior_result.lesion_energy)
    prior_increment = np.float32(
        np.float32(np.float32(np.float32(0.5) * np.float32(STANDARD13_PRIOR_COEFFICIENT))
                   * np.float32(prior_result.lesion_prior))
    )

    # PowerTable has a float lookup API but its Math.pow table is double. The
    # production Java class provides the definitive lookup semantics.
    delta = np.float32(image[1] - lesion_centroids[1])
    residual_squared = np.float32(delta * delta)
    python_power = float(toads_powertable_lookup(
        float(residual_squared), exponent, cached_exponent=exponent,
        minimum=0.0, maximum=4.0, step=float(np.float32(1.0e-6)),
    ).item())

    with tempfile.TemporaryDirectory(prefix="p41z-one-voxel-energy-") as temporary:
        work = Path(temporary)
        input_path = work / "target.bin"
        _write_input(input_path, arrays, stats, case)
        classes = work / "classes"
        classes.mkdir()
        subprocess.run([javac, "-d", str(classes), str(POWER_TABLE), str(JAVA)], check=True,
                       capture_output=True, text=True)
        java_report = json.loads(subprocess.run(
            [java, "-cp", str(classes), "P41zSingleVoxelLesionEnergyDiagnostic", str(input_path)],
            check=True, capture_output=True, text=True,
        ).stdout)

    python_report = {
        "raw_numL": float(raw), "raw_bits": _f32bits(raw),
        "smoothing_increment": float(smooth_increment), "smoothing_bits": _f32bits(smooth_increment),
        "prior_increment": float(prior_increment), "prior_bits": _f32bits(prior_increment),
        "final_energy": float(final_energy), "final_bits": _f32bits(final_energy),
        "flair_residual_squared": float(residual_squared),
        "python_powertable_double": python_power,
        "python_powertable_float32_bits": _f32bits(python_power),
    }
    terms = ("raw_numL", "smoothing_increment", "prior_increment", "final_energy")
    python_bits = {"raw_numL": "raw_bits", "smoothing_increment": "smoothing_bits",
                   "prior_increment": "prior_bits", "final_energy": "final_bits"}
    java_bits = {"raw_numL": "raw_bits", "smoothing_increment": "smoothing_bits",
                 "prior_increment": "prior_bits", "final_energy": "final_bits"}
    exact = {name: python_report[python_bits[name]] == java_report[java_bits[name]] for name in terms}
    exact["python_vs_java_powertable_raw"] = python_report["raw_bits"] == java_report["raw_bits"]
    return {
        "target_xyz": list(COORD), "class10_shape": float(shapes[k]),
        "active_lesion_class": bool(shapes[k] > 0),
        "immediately_previous_traversal_voxel_xyz": [1, 28, 28],
        "previous_class10_shape": float(arrays["shape_values"][1, 28, 28, k]),
        "python": python_report, "java_actual_powertable": java_report,
        "bitwise_match": exact, "all_terms_exact": all(exact.values()),
        "scope": "single target voxel; six fixed face neighbors; refreshed statistics recomputed from fixture; no full membership volume pass",
    }


if __name__ == "__main__":
    print(json.dumps(diagnose(), indent=2))
