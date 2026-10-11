#!/usr/bin/env python3
"""Compare one P41z iteration-6 post-statistics membership pass with Java.

The opt-in full-volume case is intentionally held fixed: iteration-6 thinning
is already complete, iteration-5 relation state is recovered from its
hash-verified checkpoint, and neither Java nor Python reruns earlier stages.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import LesionToadsThinningPass, LesionToadsThinningQueueRun
from pymacruise.lesion_next_iteration_membership import run_lesion_toads_next_gm_membership_update
from pymacruise.lesion_next_iteration_statistics import (
    LesionNextIterationStatistics,
    LesionToadsNextStatisticsRefresh,
    refresh_statistics_after_next_iteration_thinning,
)
from pymacruise.lesion_statistics import LesionToadsAuxiliaryCentroids, LesionToadsVariances
from pymacruise.lesion_standard13_volume import STANDARD13_TYPES

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests/java_reference/LesionNextIterationMembershipReferenceHarness.java"
FIXTURE = ROOT / "tests/data/synthetic/p41z_iteration6_membership_java"
STATISTICS_SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
TYPES = tuple(STANDARD13_TYPES)
OUTPUTS = (
    ("class_memberships", 13), ("class_energy", 13),
    ("lesion", 1), ("black", 1), ("outlier", 1), ("lesion_relation", 1),
    ("denominator", 1), ("lesion_energy", 1), ("black_hole_energy", 1),
    ("running_objective", 1), ("objective", 1), ("nan_count", 1),
)


def _load(path: Path, dtype: str | np.dtype | None = None) -> np.ndarray:
    image = nib.load(str(path))
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"{path} is not RAS")
    values = np.asarray(image.dataobj)
    return values.astype(dtype, copy=False) if dtype is not None else values


def _save_binary_input(path: Path, arrays: dict[str, np.ndarray], stats: LesionNextIterationStatistics) -> None:
    spatial = arrays["images"].shape[:3]
    record_type = np.dtype([
        ("image", ">f4", (2,)), ("shape", ">f4", (13,)), ("direct", "u1", (13,)),
        ("memberships", ">f4", (13,)), ("auxiliary", ">f4", (4,)), ("distances", ">i4", (3,)),
    ], align=False)
    records = np.empty(spatial, dtype=record_type)
    records["image"] = arrays["images"]
    records["shape"] = arrays["shape_values"]
    records["direct"] = arrays["direct_relations"]
    records["memberships"] = arrays["class_memberships"]
    for channel, name in enumerate(("lesion", "black", "outlier", "lesion_relation")):
        records["auxiliary"][..., channel] = arrays[name]
    records["distances"] = arrays["distance_maps"]
    with path.open("wb") as stream:
        stream.write(np.asarray(spatial, dtype=">i4").tobytes())
        state = stats.statistics
        for array in (
            state.centroids.class_centroids, state.variances.class_standard_deviation,
            state.centroids.lesion_centroids, state.variances.lesion_standard_deviation,
            state.channel_weights, np.ones(13, dtype=np.float32), np.ones(13, dtype=np.float32),
        ):
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="C"))
        stream.write(records.tobytes(order="C"))


def _stats_from_fixture() -> tuple[dict[str, np.ndarray], LesionNextIterationStatistics, dict]:
    stats_case = json.loads((FIXTURE / "case.json").read_text())
    fixture_records = stats_case.get("fixture_files", {})
    for name, metadata in fixture_records.items():
        path = FIXTURE / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["sha256"]:
            raise ValueError(f"fixture hash mismatch for {name}")
        image = nib.load(str(path))
        if nib.aff2axcodes(image.affine) != ("R", "A", "S") or list(image.shape) != metadata["shape"]:
            raise ValueError(f"fixture RAS/grid metadata mismatch for {name}")
    selected = stats_case.get("selected_java_source", {})
    if selected.get("membership_harness_sha256") != hashlib.sha256(HARNESS.read_bytes()).hexdigest():
        raise ValueError("membership harness hash does not match fixture metadata")
    inputs = stats_case["input_files"]
    images = _load(FIXTURE / inputs["images"], np.float32)
    class_memberships = _load(FIXTURE / inputs["memberships"], np.float32)
    segmentation = _load(FIXTURE / inputs["segmentation"])
    lesion = _load(FIXTURE / inputs["lesion_membership"], np.float32)
    black = _load(FIXTURE / inputs["black_hole_membership"], np.float32)
    outlier = _load(FIXTURE / inputs["outlier_membership"], np.float32)
    run = LesionToadsThinningQueueRun(segmentation, None, None, None, None, (), None)
    thinning = LesionToadsThinningPass(None, None, run)
    stats = refresh_statistics_after_next_iteration_thinning(
        thinning, images, class_memberships, lesion, black, outlier,
        stats_case["initial_class_centroids"], stats_case["initial_lesion_centroids"],
        stats_case["initial_black_hole_centroids"], TYPES,
        centroid_options={
            "mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
            "optimized_csf_all_modalities": True,
            "tissue_power_exponents": stats_case["tissue_power_exponents"],
            "lesion_power_exponents": stats_case["lesion_power_exponents"],
        },
        channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
        use_lesion_weight=False, flair_channel=1,
    )
    arrays = {
        "images": images,
        "class_memberships": class_memberships,
        "shape_values": _load(FIXTURE / "shape_values_13ch.nii.gz", np.float32),
        "direct_relations": _load(FIXTURE / "direct_relations_13ch.nii.gz", np.float32) > 0.5,
        "lesion": _load(FIXTURE / inputs["lesion_membership"], np.float32),
        "black": _load(FIXTURE / inputs["black_hole_membership"], np.float32),
        "outlier": _load(FIXTURE / inputs["outlier_membership"], np.float32),
        "lesion_relation": _load(FIXTURE / "lesion_relation.nii.gz", np.float32),
        "distance_maps": _load(FIXTURE / "post_thinning_distance_maps_3ch.nii.gz", np.int16),
    }
    if arrays["class_memberships"].shape != (*images.shape[:3], 13):
        raise ValueError("iteration-5 membership channels do not match the image grid")
    return arrays, stats, stats_case


def _parse_java(path: Path, spatial: tuple[int, int, int]) -> dict[str, np.ndarray | int]:
    # Read Java DataOutputStream's big-endian words as integers, convert those
    # words to native uint32, then reinterpret the native words as float bits.
    raw = np.fromfile(path, dtype=">u4").astype(np.uint32)
    dims = raw[:3].astype(np.uint32)
    if tuple(int(v) for v in dims) != spatial:
        raise ValueError(f"Java output dimensions mismatch: {dims}")
    total = int(np.prod(spatial))
    cursor = 3
    result: dict[str, np.ndarray | int] = {}
    for name, channels in OUTPUTS[:-2]:
        length = total * channels
        result[name] = raw[cursor:cursor + length].copy().view(np.float32).reshape((*spatial, channels) if channels > 1 else spatial)
        cursor += length
    result["objective"] = int(raw[cursor]); cursor += 1
    result["nan_count"] = int(raw[cursor]); cursor += 1
    if cursor != len(raw):
        raise ValueError("Java output has unexpected byte count")
    return result


def compare(*, java: str = "java", javac: str = "javac") -> dict:
    from compare_p41z_next_iteration_statistics_java import compare as compare_statistics

    stats_report = compare_statistics(FIXTURE, java=java, javac=javac)
    if not stats_report["all_exact"]:
        raise RuntimeError(f"upstream statistics comparison failed: {stats_report['difference_diagnostics']}")
    arrays, stats, stats_case = _stats_from_fixture()
    spatial = arrays["images"].shape[:3]
    with tempfile.TemporaryDirectory(prefix="p41z-iter6-membership-java-") as temporary:
        work = Path(temporary)
        in_path, out_path = work / "input.bin", work / "output.bin"
        _save_binary_input(in_path, arrays, stats)
        classes = work / "classes"; classes.mkdir()
        subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
        subprocess.run([java, "-Xmx4g", "-cp", str(classes), "LesionNextIterationMembershipReferenceHarness", str(in_path), str(out_path)], check=True, capture_output=True, text=True)
        reference = _parse_java(out_path, spatial)

    result = run_lesion_toads_next_gm_membership_update(
        stats, arrays["images"], arrays["shape_values"], np.ones(13, np.float32), np.ones(13, np.float32),
        arrays["direct_relations"], arrays["class_memberships"], arrays["lesion"], arrays["black"],
        arrays["outlier"], arrays["lesion_relation"], arrays["distance_maps"][..., 0],
        arrays["distance_maps"][..., 1], arrays["distance_maps"][..., 2],
        previous_energy=0.0, outlier_energy=0.1, lesion_power_exponent=1.2,
        spread=0.5, source_zero_variance_nan=True, intensity_fields=None,
    )
    actual = result.membership
    pairs = {
        "class_memberships": actual.class_memberships,
        "class_energy": actual.class_energies,
        "lesion": actual.lesion_membership,
        "black": actual.black_hole_membership,
        "outlier": actual.outlier_membership,
        "lesion_relation": actual.lesion_relation,
        "denominator": actual.denominator,
        "lesion_energy": actual.lesion_energy,
        "black_hole_energy": actual.black_hole_energy,
        "running_objective": actual.running_objective,
    }
    exact, nan_counts, differences = {}, {"java": int(reference["nan_count"]), "python": 0}, {}
    for name, values in pairs.items():
        expected = np.asarray(reference[name], dtype=np.float32)
        values = np.asarray(values, dtype=np.float32)
        nan_counts["python"] += int(np.isnan(values).sum())
        diff = values.view(np.uint32) != expected.view(np.uint32)
        exact[name] = not bool(np.any(diff))
        if np.any(diff):
            finite = np.isfinite(values) & np.isfinite(expected)
            delta = np.abs(values.astype(np.float64) - expected.astype(np.float64))
            differences[name] = {
                "different": int(diff.sum()),
                "first_indices": np.argwhere(diff)[:10].tolist(),
                "max_finite_abs_error": float(delta[finite].max()) if finite.any() else None,
            }
    expected_objective = int(reference["objective"])
    objective_exact = int(np.asarray(actual.objective, np.float32).view(np.uint32)) == expected_objective
    exact["objective"] = objective_exact
    nan_counts["python"] += int(np.isnan(actual.objective))
    return {
        "spatial_shape_xyz": list(spatial),
        "profile": {"identity_weights": True, "outlier_energy": 0.1, "lesion_power_exponent": 1.2,
                    "use_lesion_weight": True, "gm_flag": True, "fields": None},
        "statistics_exact_first": stats_report["all_exact"],
        "compared_stage_order": [*stats.stage_order, "computeMemberships:GMflag=true"],
        "energy_ratio_checked": False,
        "nan_counts": nan_counts,
        "exact_float32_bits": exact,
        "all_exact": all(exact.values()) and nan_counts["java"] == nan_counts["python"],
        "difference_diagnostics": differences,
        "java_objective_bits": expected_objective,
        "python_objective_bits": int(np.asarray(actual.objective, np.float32).view(np.uint32)),
        "source_sha256": STATISTICS_SOURCE_SHA,
        "wrapper_sha256": WRAPPER_SHA,
        "java_harness_sha256": hashlib.sha256(HARNESS.read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--java", default="java")
    parser.add_argument("--javac", default="javac")
    args = parser.parse_args()
    report = compare(java=args.java, javac=args.javac)
    print(json.dumps(report, indent=2))
    if not report["all_exact"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
