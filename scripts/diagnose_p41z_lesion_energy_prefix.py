#!/usr/bin/env python3
"""Compare P41z lesion-energy terms after the traversal prefix to [1,28,29]."""
from __future__ import annotations

import json
import struct
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from compare_p41z_iteration6_membership_java import (
    FIXTURE, HARNESS, _save_binary_input, _stats_from_fixture,
)
from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies
from pymacruise.lesion_membership import (
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
)
from pymacruise.lesion_similarity import compute_lesion_toads_similarity
from pymacruise.lesion_standard13_voxel import (
    STANDARD13_CLASS_SMOOTHING, STANDARD13_NAMES, STANDARD13_PRIOR_COEFFICIENT,
    STANDARD13_PRIOR_SCALE, STANDARD13_TYPES, compute_lesion_toads_standard13_voxel,
)
from pymacruise.toads_math import toads_powertable_lookup

ROOT = Path(__file__).resolve().parents[1]
JAVA_PREFIX = ROOT / "tests/java_reference/P41zTraversalPrefixCapture.java"
JAVA_TERMS = ROOT / "tests/java_reference/P41zSingleVoxelLesionEnergyDiagnostic.java"
TOADS = ROOT.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic"
POWER_TABLE = TOADS / "utilities/PowerTable.java"
COORD = (1, 28, 29)
FACES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
CROP = (slice(0, 3), slice(0, 30), slice(0, 31))


def _bits(value) -> str:
    return struct.pack(">f", np.float32(value)).hex()


def _face_state(arrays: dict[str, np.ndarray], point: tuple[int, int, int]):
    x, y, z = point
    classes = np.asarray([arrays["class_memberships"][x+dx, y+dy, z+dz] for dx,dy,dz in FACES], dtype=np.float32)
    lesion = np.asarray([arrays["lesion"][x+dx, y+dy, z+dz] for dx,dy,dz in FACES], dtype=np.float32)
    black = np.asarray([arrays["black"][x+dx, y+dy, z+dz] for dx,dy,dz in FACES], dtype=np.float32)
    return classes, lesion, black


def _python_terms(arrays, stats, case, face_state):
    x, y, z = COORD
    settings = case["membership_settings"]
    k = 10
    image = arrays["images"][x, y, z]
    shapes = arrays["shape_values"][x, y, z]
    face_classes, face_lesion, face_black = face_state
    centers = stats.statistics.centroids.class_centroids
    deviations = stats.statistics.variances.class_standard_deviation
    lesion_center = stats.statistics.centroids.lesion_centroids
    lesion_deviation = stats.statistics.variances.lesion_standard_deviation
    weights = stats.statistics.channel_weights
    exponent = float(settings["lesion_power_exponent"])
    raw_result = compute_lesion_toads_unfielded_data_energies(
        image, centers[:,k], lesion_center, np.zeros(2, np.float32), weights,
        object_type="lesion", flair_channel=1,
        use_lesion_weight=bool(settings["use_lesion_weight"]),
        lesion_power_exponent=exponent, cached_lesion_power_exponent=exponent,
    )
    raw = np.float32(raw_result.lesion_energy)
    smooth_result = update_lesion_toads_smoothing_energies(
        face_classes, face_lesion, face_black, class_index=k, class_type="LESION",
        smoothing_prior=float(settings["smoothing_priors"][k]),
        smoothing=float(STANDARD13_CLASS_SMOOTHING), class_energy=0.0,
        lesion_energy=float(raw), black_hole_energy=0.0,
        class_membership_data_term=0.0, lesion_membership_data_term=float(raw),
        black_hole_membership_data_term=0.0, outmask_membership_data_term=0.0,
    )
    smooth_sum = np.float32(smooth_result.other_class_neighbor_squares
                            + smooth_result.current_class_neighbor_squares)
    smooth_sum = np.float32(smooth_sum + smooth_result.black_hole_neighbor_squares)
    smooth_increment = np.float32(
        np.float32(np.float32(np.float32(0.5)*np.float32(settings["smoothing_priors"][k]))
                   *np.float32(STANDARD13_CLASS_SMOOTHING))*smooth_sum
    )
    similarities = compute_lesion_toads_similarity(
        centers, deviations, lesion_center, lesion_deviation, weights,
        STANDARD13_TYPES, STANDARD13_NAMES, prior_scale=float(STANDARD13_PRIOR_SCALE),
        t1_channel=0, black_hole_similarity=np.zeros(13,np.float32),
        finite_zero_variance=not bool(settings["source_zero_variance_nan"]),
    )
    prior_result = update_lesion_toads_prior_energies(
        np.asarray(settings["shape_priors"],np.float32), similarities.class_similarity[k],
        similarities.lesion_similarity, similarities.black_hole_similarity, shapes,
        class_index=k, class_type="LESION", prior_coefficient=float(STANDARD13_PRIOR_COEFFICIENT),
        class_membership_data_term=0.0,
        lesion_membership_data_term=float(smooth_result.lesion_membership_data_term),
        black_hole_membership_data_term=0.0, outmask_membership_data_term=0.0,
        class_energy=0.0, lesion_energy=float(smooth_result.lesion_energy), black_hole_energy=0.0,
        has_lesion=True, lesion_class_index=10, lesion_class_type="LESION",
    )
    prior_increment = np.float32(
        np.float32(np.float32(np.float32(0.5)*np.float32(STANDARD13_PRIOR_COEFFICIENT))
                   *np.float32(prior_result.lesion_prior))
    )
    final_energy = np.float32(prior_result.lesion_energy)
    residual = np.float32(image[1]-lesion_center[1])
    residual_squared = np.float32(residual*residual)
    powertable = float(toads_powertable_lookup(
        float(residual_squared), exponent, cached_exponent=exponent,
        minimum=0.0, maximum=4.0, step=float(np.float32(1e-6)),
    ).item())
    return {
        "raw_bits": _bits(raw), "smoothing_bits": _bits(smooth_increment),
        "prior_bits": _bits(prior_increment), "final_bits": _bits(final_energy),
        "raw": float(raw), "smoothing_increment": float(smooth_increment),
        "prior_increment": float(prior_increment), "final_energy": float(final_energy),
        "python_powertable_double": powertable,
    }


def diagnose(java="java", javac="javac"):
    full, stats, case = _stats_from_fixture()
    arrays = {name: value[CROP].copy() for name, value in full.items()}
    # Keep the Java harness input deliberately small: x=0:3, y=0:30, z=0:31.
    with tempfile.TemporaryDirectory(prefix="p41z-traversal-prefix-") as temporary:
        work = Path(temporary)
        input_path, capture_path = work/"crop.bin", work/"target.bin"
        _save_binary_input(input_path, arrays, stats)
        classes = work/"classes"; classes.mkdir()
        subprocess.run([javac,"-d",str(classes),str(POWER_TABLE),str(HARNESS),
                        str(JAVA_TERMS),str(JAVA_PREFIX)],check=True,capture_output=True,text=True)
        java_prefix = json.loads(subprocess.run(
            [java,"-cp",str(classes),"P41zTraversalPrefixCapture",str(input_path),str(capture_path)],
            check=True,capture_output=True,text=True,
        ).stdout)
        java_terms = json.loads(subprocess.run(
            [java,"-cp",str(classes),"P41zSingleVoxelLesionEnergyDiagnostic",str(capture_path)],
            check=True,capture_output=True,text=True,
        ).stdout)

    settings = case["membership_settings"]
    shapes = arrays["shape_values"]
    current = arrays["class_memberships"]
    lesion = arrays["lesion"].copy(); black = arrays["black"].copy()
    outlier = arrays["outlier"].copy(); lesion_relation = arrays["lesion_relation"].copy()
    direct = arrays["direct_relations"]
    distances = arrays["distance_maps"]
    state_energy = np.zeros(13,np.float32); lesion_energy_state = np.float32(0.0)
    similarity = compute_lesion_toads_similarity(
        stats.statistics.centroids.class_centroids, stats.statistics.variances.class_standard_deviation,
        stats.statistics.centroids.lesion_centroids, stats.statistics.variances.lesion_standard_deviation,
        stats.statistics.channel_weights, STANDARD13_TYPES, STANDARD13_NAMES,
        prior_scale=float(STANDARD13_PRIOR_SCALE), t1_channel=0,
        black_hole_similarity=np.zeros(13,np.float32),
        finite_zero_variance=not bool(settings["source_zero_variance_nan"]),
    )
    python_face = None
    target_result = None
    tx,ty,tz = COORD
    for x in range(1,2):
        for y in range(1,ty+1):
            for z in range(1,(tz+1 if y==ty else arrays["images"].shape[2]-1)):
                if (x,y,z)==COORD:
                    python_face = _face_state({"class_memberships":current,"lesion":lesion,"black":black},COORD)
                faces = _face_state({"class_memberships":current,"lesion":lesion,"black":black},(x,y,z))
                result = compute_lesion_toads_standard13_voxel(
                    arrays["images"][x,y,z], stats.statistics.centroids.class_centroids,
                    stats.statistics.variances.class_standard_deviation,
                    stats.statistics.centroids.lesion_centroids,
                    stats.statistics.variances.lesion_standard_deviation,
                    stats.statistics.channel_weights, shapes[x,y,z],
                    np.asarray(settings["shape_priors"],np.float32),
                    np.asarray(settings["smoothing_priors"],np.float32),direct[x,y,z],
                    current[x,y,z],lesion[x,y,z],black[x,y,z],outlier[x,y,z],lesion_relation[x,y,z],
                    *faces, outlier_energy=float(settings["outlier_energy"]),
                    lesion_power_exponent=float(settings["lesion_power_exponent"]),
                    ventricle_distance=int(distances[x,y,z,0]),
                    gray_matter_distance=int(distances[x,y,z,1]),
                    brainstem_distance=int(distances[x,y,z,2]),gm_flag=True,
                    spread=float(settings["spread"]),
                    source_zero_variance_nan=bool(settings["source_zero_variance_nan"]),
                    precomputed_similarities=similarity,
                    initial_class_energies=state_energy,
                    initial_lesion_energy=float(lesion_energy_state),
                )
                current[x,y,z]=result.class_memberships
                lesion[x,y,z]=result.lesion_membership; black[x,y,z]=result.black_hole_membership
                outlier[x,y,z]=result.outlier_membership; lesion_relation[x,y,z]=result.lesion_relation
                state_energy=result.class_energies.copy(); lesion_energy_state=result.lesion_energy
                if (x,y,z)==COORD:
                    target_result=result
                    break
            if target_result is not None: break
        if target_result is not None: break
    if python_face is None or target_result is None: raise RuntimeError("target traversal was not reached")

    python_terms = _python_terms(arrays,stats,case,python_face)
    face_equal = (
        [_bits(v) for v in python_face[0].ravel()] == java_prefix["face_class_bits"]
        and [_bits(v) for v in python_face[1]] == java_prefix["face_lesion_bits"]
        and [_bits(v) for v in python_face[2]] == java_prefix["face_black_bits"]
    )
    term_match = {
        "raw_numL": python_terms["raw_bits"] == java_terms["raw_bits"],
        "smoothing_increment": python_terms["smoothing_bits"] == java_terms["smoothing_bits"],
        "prior_increment": python_terms["prior_bits"] == java_terms["prior_bits"],
        "final_energy": python_terms["final_bits"] == java_terms["final_bits"],
        "target_voxel_energy_after_prefix": _bits(target_result.lesion_energy) == java_terms["final_bits"],
    }
    return {
        "crop_xyz": [3,30,31], "target_xyz": list(COORD),
        "updated_prefix_voxels": 1*(ty-1)*(arrays["images"].shape[2]-2)+(tz-1),
        "python_face_class_bits": [_bits(v) for v in python_face[0].ravel()],
        "java_face_class_bits": java_prefix["face_class_bits"],
        "python_face_lesion_bits": [_bits(v) for v in python_face[1]],
        "java_face_lesion_bits": java_prefix["face_lesion_bits"],
        "python_face_black_bits": [_bits(v) for v in python_face[2]],
        "java_face_black_bits": java_prefix["face_black_bits"],
        "six_face_state_exact": face_equal,
        "python_terms": python_terms,
        "java_terms": java_terms,
        "python_target_final_energy_bits": _bits(target_result.lesion_energy),
        "java_source_decomposition_final_energy_bits": java_terms["final_bits"],
        "existing_java_oracle_map_energy_bits": java_prefix["target_final_lesion_energy_bits"],
        "existing_java_oracle_map_energy_matches_source_decomposition":
            java_prefix["target_final_lesion_energy_bits"] == java_terms["final_bits"],
        "term_bitwise_matches": term_match,
        "all_exact": face_equal and all(term_match.values()),
        "scope": "3x30x31 cropped fixture; in-place traversal stops immediately before and at [1,28,29]",
    }


if __name__ == "__main__":
    print(json.dumps(diagnose(),indent=2))
