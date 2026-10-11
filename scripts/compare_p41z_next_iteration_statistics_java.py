#!/usr/bin/env python3
"""Compare P41z next-loop statistics against a selected Java source excerpt.

Input directory contract (all images are RAS NIfTI):
  case.json: channels, class_types, initial centroids, auxiliary centroids,
             class/lesion power exponents, channel norm type, and outlier
             distance.
  images.nii.gz: (x,y,z,2) float32 [T1_MPRAGE, FLAIR]
  memberships.nii.gz: (x,y,z,13) float32
  segmentation.nii.gz: (x,y,z) integer
  lesion_membership.nii.gz, black_hole_membership.nii.gz,
  outlier_membership.nii.gz: (x,y,z) float32

The script serializes those NIfTIs into the Java harness's documented binary
contract, runs the source excerpt, invokes the Python production seam, and
reports exact float-bit parity plus per-side NaN counts. The harness uses a
fixed selected-source category profile, Identity channel weighting, and no
correction fields.
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

from pymacruise.lesion_growth import (
    LesionToadsThinningPass,
    LesionToadsThinningQueueRun,
)
from pymacruise.lesion_next_iteration_statistics import (
    refresh_statistics_after_next_iteration_thinning,
)


ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests/java_reference/LesionNextIterationStatisticsReferenceHarness.java"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
TYPES = ("mask", "optimized_csf", "optimized_csf", "object", "object", "object", "object", "object", "object", "object", "lesion", "object", "object")
TYPE_CODES = {"mask": 0, "optimized_csf": 1, "object": 2, "lesion": 3}


def _load(path: Path, dtype: np.dtype | None = None) -> np.ndarray:
    image = nib.load(str(path))
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"{path.name} must be RAS; got {nib.aff2axcodes(image.affine)}")
    values = np.asarray(image.dataobj)
    return values.astype(dtype, copy=False) if dtype is not None else values


def _write_array(stream, array: np.ndarray, dtype: str) -> None:
    stream.write(np.asarray(array, dtype=dtype).tobytes(order="C"))


def _read_bits(value: list, shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(value, dtype=np.uint32).reshape(shape).view(np.float32)


def compare(input_dir: Path, *, java: str = "java", javac: str = "javac") -> dict:
    case = json.loads((input_dir / "case.json").read_text(encoding="utf-8"))
    input_files = case.get("input_files", {})
    paths = {
        "images": input_files.get("images", "images.nii.gz"),
        "memberships": input_files.get("memberships", "memberships.nii.gz"),
        "segmentation": input_files.get("segmentation", "segmentation.nii.gz"),
        "lesion_membership": input_files.get("lesion_membership", "lesion_membership.nii.gz"),
        "black_hole_membership": input_files.get("black_hole_membership", "black_hole_membership.nii.gz"),
        "outlier_membership": input_files.get("outlier_membership", "outlier_membership.nii.gz"),
    }
    fixture_files = case.get("fixture_files", {})
    for name, metadata in fixture_files.items():
        path = input_dir / name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != metadata["sha256"]:
            raise ValueError(f"fixture hash mismatch for {name}")
        image = nib.load(str(path))
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"fixture image is not RAS: {name}")
        if list(image.shape) != metadata["shape"] or str(image.get_data_dtype()) != metadata["dtype"]:
            raise ValueError(f"fixture header metadata mismatch for {name}")
    images = _load(input_dir / paths["images"], np.dtype(np.float32))
    memberships = _load(input_dir / paths["memberships"], np.dtype(np.float32))
    segmentation = _load(input_dir / paths["segmentation"], np.dtype(np.int32))
    lesion = _load(input_dir / paths["lesion_membership"], np.dtype(np.float32))
    black = _load(input_dir / paths["black_hole_membership"], np.dtype(np.float32))
    outlier = _load(input_dir / paths["outlier_membership"], np.dtype(np.float32))
    spatial = segmentation.shape
    classes = len(TYPES)
    if images.shape != (*spatial, 2) or memberships.shape != (*spatial, classes):
        raise ValueError("expected images (x,y,z,2) and memberships (x,y,z,13)")
    if any(array.shape != spatial for array in (lesion, black, outlier)):
        raise ValueError("all auxiliary maps must match segmentation grid")
    if tuple(case["object_types"]) != TYPES:
        raise ValueError("case object_types must match the pinned 13-class P41z profile")
    if case.get("norm_type") != "Identity" or float(case.get("outlier_distance", -1)) != 0.1:
        raise ValueError("P41z comparison requires Identity normalization and outlier_distance 0.1")
    if not np.array_equal(np.asarray(case["tissue_power_exponents"], dtype=np.float32), np.ones((2, 13), dtype=np.float32)):
        raise ValueError("P41z comparison requires unit tissue power exponents")
    if not np.array_equal(np.asarray(case["lesion_power_exponents"], dtype=np.float32), np.asarray([0.4, 1.2], dtype=np.float32)):
        raise ValueError("P41z comparison requires lesion powers (0.4, 1.2)")
    if case.get("optimized_csf_all_modalities") is not True:
        raise ValueError("P41z VandyCRUISE profile requires optimized CSF PowerTable on every modality")
    selected = case.get("selected_java_source", {})
    if selected and (
        selected.get("source_sha256") != SOURCE_SHA
        or selected.get("wrapper_sha256") != WRAPPER_SHA
    ):
        raise ValueError("fixture selected Java source hashes do not match this comparator")
    if selected.get("harness_sha256"):
        harness_hash = hashlib.sha256(HARNESS.read_bytes()).hexdigest()
        if selected["harness_sha256"] != harness_hash:
            raise ValueError("fixture Java harness hash does not match current source")
    nx, ny, nz = spatial
    nc = 2
    with tempfile.TemporaryDirectory(prefix="p41z-next-stats-java-") as td:
        work = Path(td)
        input_path, output_path = work / "input.bin", work / "reference.json"
        with input_path.open("wb") as stream:
            stream.write(np.asarray([nx, ny, nz, nc, classes], dtype=">i4").tobytes())
            stream.write(np.asarray([TYPE_CODES[kind] for kind in TYPES], dtype=np.uint8).tobytes())
            _write_array(stream, images, ">f4")
            _write_array(stream, memberships, ">f4")
            _write_array(stream, segmentation, ">i4")
            for values in (lesion, black, outlier):
                _write_array(stream, values, ">f4")
            _write_array(stream, np.asarray(case["initial_class_centroids"], dtype=np.float32), ">f4")
            _write_array(stream, np.asarray(case["initial_lesion_centroids"], dtype=np.float32), ">f4")
            _write_array(stream, np.asarray(case["initial_black_hole_centroids"], dtype=np.float32), ">f4")
            _write_array(stream, np.asarray(case["tissue_power_exponents"], dtype=np.float32), ">f4")
            _write_array(stream, np.asarray(case["lesion_power_exponents"], dtype=np.float32), ">f4")
        classes_dir = work / "classes"
        classes_dir.mkdir()
        subprocess.run([javac, "-d", str(classes_dir), str(HARNESS)], check=True, capture_output=True, text=True)
        subprocess.run([java, "-cp", str(classes_dir), "LesionNextIterationStatisticsReferenceHarness", str(input_path), str(output_path)], check=True, capture_output=True, text=True)
        reference = json.loads(output_path.read_text(encoding="utf-8"))

    # The seam consumes a thinning-pass object but reads only its run's
    # segmentation. Unrelated queue state is intentionally absent here.
    run = LesionToadsThinningQueueRun(segmentation, None, None, None, None, (), None)
    thinning = LesionToadsThinningPass(None, None, run)
    result = refresh_statistics_after_next_iteration_thinning(
        thinning,
        images,
        memberships,
        lesion,
        black,
        outlier,
        case["initial_class_centroids"],
        case["initial_lesion_centroids"],
        case["initial_black_hole_centroids"],
        TYPES,
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "optimized_csf_all_modalities": True,
            "tissue_power_exponents": case["tissue_power_exponents"],
            "lesion_power_exponents": case["lesion_power_exponents"],
        },
        channel_weight_options={"outlier_distance": case["outlier_distance"], "norm_type": case["norm_type"]},
        use_lesion_weight=False,
        flair_channel=1,
    )
    stats = result.statistics
    actual = {
        "class_centroid_bits": np.asarray(stats.centroids.class_centroids, dtype=np.float32),
        "lesion_centroid_bits": np.asarray(stats.centroids.lesion_centroids, dtype=np.float32),
        "black_hole_centroid_bits": np.asarray(stats.centroids.black_hole_centroids, dtype=np.float32),
        "class_std_bits": np.asarray(stats.variances.class_standard_deviation, dtype=np.float32),
        "lesion_std_bits": np.asarray(stats.variances.lesion_standard_deviation, dtype=np.float32),
        "black_hole_std_bits": np.asarray(stats.variances.black_hole_standard_deviation, dtype=np.float32),
        "channel_weight_bits": np.asarray(stats.channel_weights, dtype=np.float32),
    }
    expected_shapes = {
        "class_centroid_bits": (nc, classes), "lesion_centroid_bits": (nc,),
        "black_hole_centroid_bits": (nc,), "class_std_bits": (nc, classes),
        "lesion_std_bits": (nc,), "black_hole_std_bits": (nc,),
        "channel_weight_bits": (nc,),
    }
    exact = {}
    diagnostics = {}
    counts = {"java": int(reference["nan_count"]), "python": 0}
    for key, values in actual.items():
        values = values.reshape(expected_shapes[key])
        java_values = _read_bits(reference[key], expected_shapes[key])
        counts["python"] += int(np.isnan(values).sum())
        different = values.view(np.uint32) != java_values.view(np.uint32)
        exact[key] = not bool(np.any(different))
        if np.any(different):
            finite = np.isfinite(values) & np.isfinite(java_values)
            delta = np.abs(values.astype(np.float64) - java_values.astype(np.float64))
            diagnostics[key] = {
                "differing_values": int(np.count_nonzero(different)),
                "max_abs_error_finite": float(np.max(delta[finite])) if np.any(finite) else None,
                "first_indices": np.argwhere(different)[:12].tolist(),
            }
        actual[key] = values
    return {
        "source_sha256": SOURCE_SHA,
        "wrapper_sha256": WRAPPER_SHA,
        "spatial_shape_xyz": list(spatial),
        "stage_order": list(result.stage_order),
        "norm_type": case["norm_type"],
        "outlier_distance": case["outlier_distance"],
        "optimized_csf_all_modalities": case["optimized_csf_all_modalities"],
        "nan_counts": counts,
        "exact_float32_bits": exact,
        "difference_diagnostics": diagnostics,
        "all_exact": all(exact.values()),
        "java_statistic_bits": {
            key: reference[key]
            for key in actual
        },
        "java_nan_channel_weight_bits": reference["channel_weight_bits"],
        "python_channel_weight_bits": actual["channel_weight_bits"].view(np.uint32).astype(np.uint64).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--java", default="java")
    parser.add_argument("--javac", default="javac")
    args = parser.parse_args()
    print(json.dumps(compare(args.input_dir, java=args.java, javac=args.javac), indent=2))


if __name__ == "__main__":
    main()
