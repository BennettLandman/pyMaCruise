#!/usr/bin/env python3
"""Run selected JIST CubicVolumeCropper on a saved human base-prep attempt.

This calls the compiled upstream JIST cropper class directly. The Java side
derives one crop from Python's saved MaSeg-derived Atlas_GM, then applies that
same CropParameters instance to the selected target arrays. It does not rerun
transfer, cortical separation, fusion, or any later CRUISE stages.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any
import zipfile

import nibabel as nib
import numpy as np


REPO = Path(__file__).resolve().parents[1]
MIPAV = REPO.parent / "mipav-src"
JIST = REPO.parent / "jist-scm-2026-09-26"
HARNESS = REPO / "tests" / "java_reference" / "BaseTissueCropperReferenceHarness.java"
WRAPPER_SOURCE = (
    REPO.parent / "toads-cruise-scm-2026-09-26" / "toads-cruise" / "src"
    / "edu" / "vanderbilt" / "masi" / "plugins" / "CRUISE" / "utilities"
    / "CropBasedOnGM.java"
)
RUNTIME_JAR = REPO / "build" / "java-reference" / "jist-runtime-without-selected-conflicts.jar"
BUILD = REPO / "build" / "java-reference" / "human-base-cropper"

JAVA_SUFFIX = {True: ".i32be", False: ".f32be"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jar_entry_sha256(jar_path: Path, entry: str) -> str:
    with zipfile.ZipFile(jar_path) as archive:
        return hashlib.sha256(archive.read(entry)).hexdigest()


def classpath() -> str:
    entries = [
        str(BUILD / "classes"),
        str(REPO / "build" / "java-reference" / "vandycruise-core-classes"),
        str(REPO / "build" / "java-reference" / "toads-orientation-mapping-classes"),
        str(REPO / "build" / "java-reference" / "clean-source-compile"),
        str(RUNTIME_JAR),
    ]
    for library_root in (
        MIPAV / "src" / "lib",
        MIPAV / "src" / "libnative",
        JIST / "JIST" / "lib",
    ):
        if library_root.exists():
            entries.extend(str(path) for path in sorted(library_root.rglob("*.jar")))
    return os.pathsep.join(entries)


def compile_harness() -> str:
    if not HARNESS.is_file():
        raise FileNotFoundError(f"Missing Java reference harness: {HARNESS}")
    if not RUNTIME_JAR.is_file():
        raise FileNotFoundError(f"Missing compiled selected JIST runtime: {RUNTIME_JAR}")
    BUILD.mkdir(parents=True, exist_ok=True)
    classes = BUILD / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    cp = classpath()
    subprocess.run(
        ["javac", "-classpath", cp, "-d", str(classes), str(HARNESS)],
        check=True,
    )
    return cp


def resolve_attempt(run_root: Path, run_data: dict[str, Any],
                    attempt: int | None) -> tuple[int, Path]:
    stage_root = run_root / "stages" / "base-tissue-preparation"
    stage_record = run_data.get("stages", {}).get("base-tissue-preparation", {})
    if attempt is None:
        if stage_record.get("status") != "complete" or not isinstance(stage_record.get("attempt"), int):
            raise RuntimeError(
                "No completed base-tissue-preparation attempt is recorded in run.json; "
                "wait until the stage completes or choose a completed run root"
            )
        attempt = int(stage_record["attempt"])
        directory = stage_root / f"attempt-{attempt:02d}"
    else:
        directory = stage_root / f"attempt-{attempt:02d}"
        if not directory.is_dir():
            raise FileNotFoundError(f"No selected preparation attempt: {directory}")
        if stage_record.get("attempt") != attempt or stage_record.get("status") != "complete":
            raise RuntimeError(
                f"run.json does not record base-tissue-preparation attempt {attempt} as complete"
            )
    if not directory.is_dir():
        raise FileNotFoundError(f"Selected completed attempt directory is missing: {directory}")
    if not (directory / "provenance.json").is_file():
        raise FileNotFoundError(f"Attempt has no provenance.json: {directory}")
    return int(attempt), directory


def resolved_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def record_input(path: Path) -> dict[str, Any]:
    return {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}


def verify_manifest_hash(path: Path, manifest_records: list[dict[str, Any]]) -> None:
    resolved = path.resolve()
    for record in manifest_records:
        recorded_path = Path(record.get("path", "")).expanduser()
        if recorded_path.resolve() == resolved:
            actual = sha256(path)
            if actual != record.get("sha256"):
                raise ValueError(f"Workflow manifest hash mismatch for {path}")
            return
    raise ValueError(f"Workflow run manifest does not record selected stage output: {path}")


def nifti(path: Path) -> nib.spatialimages.SpatialImage:
    image = nib.load(str(path))
    if len(image.shape) != 3:
        raise ValueError(f"Expected 3D image, got {image.shape}: {path}")
    return image


def assert_same_grid(reference: nib.spatialimages.SpatialImage,
                     candidate: nib.spatialimages.SpatialImage, name: str) -> None:
    if reference.shape != candidate.shape:
        raise ValueError(f"{name} shape {candidate.shape} differs from crop-mask shape {reference.shape}")
    if not np.allclose(reference.affine, candidate.affine, rtol=0, atol=1e-5):
        raise ValueError(f"{name} affine differs from crop-mask affine")


def pack_raw(image: nib.spatialimages.SpatialImage, path: Path, integer: bool) -> None:
    data = np.asanyarray(image.dataobj)
    if not np.isfinite(data).all():
        raise ValueError(f"Non-finite input voxels: {path}")
    if integer:
        if not np.equal(data, np.floor(data)).all():
            raise ValueError(f"Non-integer input voxels for integer Java adapter: {path}")
        packed = np.asarray(data, dtype=">i4", order="C")
    else:
        packed = np.asarray(data, dtype=">f4", order="C")
    packed.tofile(path)


def java_image_data(path: Path, integer: bool) -> np.ndarray:
    dtype = ">i4" if integer else ">f4"
    values = np.fromfile(path, dtype=dtype)
    return values


def compare_arrays(java: np.ndarray, python: np.ndarray, name: str,
                   expected_shape: tuple[int, int, int]) -> dict[str, Any]:
    expected = np.asarray(python)
    if java.size != expected.size or expected.shape != expected_shape:
        return {
            "name": name,
            "java_voxels": int(java.size),
            "python_voxels": int(expected.size),
            "different_voxels": None,
            "python_shape": list(expected.shape),
            "expected_shape": list(expected_shape),
            "shape_match": bool(java.size == expected.size and expected.shape == expected_shape),
        }
    java64 = java.astype(np.float64, copy=False)
    expected64 = expected.astype(np.float64, copy=False).reshape(-1)
    differences = java64 != expected64
    absolute = np.abs(java64 - expected64)
    return {
        "name": name,
        "shape_match": True,
        "python_shape": list(expected.shape),
        "expected_shape": list(expected_shape),
        "voxel_count": int(expected.size),
        "different_voxels": int(np.count_nonzero(differences)),
        "max_absolute_error": float(absolute.max(initial=0.0)),
        "mean_absolute_error": float(absolute.mean()) if absolute.size else 0.0,
        "java_numeric_type": str(java.dtype),
        "python_numeric_type": str(expected.dtype),
        "exact_values_match": bool(np.array_equal(java64, expected64)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True,
                        help="restartable base workflow root containing run.json")
    parser.add_argument("--attempt", type=int,
                        help="completed base-tissue-preparation attempt number (default: completed attempt recorded in run.json)")
    parser.add_argument("--output-dir", type=Path,
                        help="report directory (default: build/java-reference/human-base-cropper/<run>-attempt-NN)")
    parser.add_argument("--java-heap", default="1g",
                        help="Java heap cap (default: 1g)")
    args = parser.parse_args()

    run_root = args.run_root.expanduser().resolve()
    run_manifest = run_root / "run.json"
    if not run_manifest.is_file():
        raise FileNotFoundError(f"Missing workflow run manifest: {run_manifest}")
    run_data = json.loads(run_manifest.read_text())
    attempt_number, attempt_dir = resolve_attempt(run_root, run_data, args.attempt)
    stage_record = run_data["stages"]["base-tissue-preparation"]
    provenance_path = attempt_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    input_paths = provenance.get("inputs", {})
    maseg_path = resolved_path(input_paths["maseg"], run_root)
    if not maseg_path.is_file():
        raise FileNotFoundError(f"MaSeg input from preparation provenance is missing: {maseg_path}")
    run_input_record = run_data.get("inputs", {}).get("maseg", {})
    if run_input_record.get("sha256") != sha256(maseg_path):
        raise ValueError("Selected MaSeg does not match the hash recorded by run.json")

    mask_path = attempt_dir / "atlas_gm.nii.gz"
    if not mask_path.is_file():
        raise FileNotFoundError(f"Missing bounding mask in selected attempt: {mask_path}")
    gm_mask = nifti(mask_path)
    maseg = nifti(maseg_path)
    assert_same_grid(gm_mask, maseg, "MaSeg")

    targets: dict[str, tuple[Path, bool, str]] = {
        "atlas_gm": (mask_path, True, "atlas_gm_cropped"),
        "atlas_wm": (attempt_dir / "atlas_wm.nii.gz", True, "atlas_wm_cropped"),
        "atlas_brain_tissue_mask": (
            attempt_dir / "atlas_brain_tissue_mask.nii.gz", True,
            "atlas_brain_tissue_mask_cropped",
        ),
        "toads_gm": (attempt_dir / "toads_gm_transferred.nii.gz", False, "toads_gm_cropped"),
        "toads_wmfill": (
            attempt_dir / "toads_wmfill_transferred.nii.gz", False,
            "toads_wmfill_cropped",
        ),
        "toads_mask": (
            attempt_dir / "toads_mask_transferred.nii.gz", False,
            "toads_mask_cropped",
        ),
        "maseg": (maseg_path, True, ""),
    }
    for name, (path, _integer, _output_name) in targets.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing required source-transfer input {name}: {path}")
        assert_same_grid(gm_mask, nifti(path), name)
        if name != "maseg":
            verify_manifest_hash(path, stage_record.get("outputs", []))

    parameters = provenance.get("parameters", {}).get("crop", {})
    crop_edge = int(parameters.get("edge_width_voxels", run_data.get("options", {}).get("crop_edge_width", 5)))
    manifest_bounds = tuple(tuple(int(v) for v in axis) for axis in parameters["bounds_inclusive"])
    shape = tuple(int(v) for v in gm_mask.shape)
    python_affine = np.asarray(parameters["original_affine"], dtype=np.float64)
    if tuple(int(v) for v in parameters.get("original_shape", ())) != shape:
        raise ValueError("Saved crop provenance shape differs from Atlas_GM input")
    if not np.allclose(gm_mask.affine, python_affine, rtol=0.0, atol=1e-5):
        raise ValueError("Saved crop provenance affine differs from Atlas_GM input")

    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else BUILD / f"{run_root.name}-attempt-{attempt_number:02d}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    java_output = output_dir / "java-output"

    with tempfile.TemporaryDirectory(prefix="pymacruise-crop-", dir=output_dir) as temp_name:
        raw_input = Path(temp_name) / "input"
        raw_input.mkdir()
        java_output = Path(temp_name) / "java-output"
        java_home = Path(temp_name) / "java-home"
        java_home.mkdir()
        for name, (path, integer, _output_name) in targets.items():
            pack_raw(nifti(path), raw_input / f"{name}{JAVA_SUFFIX[integer]}", integer)

        cp = compile_harness()
        classpath_value = cp
        java_log = output_dir / "java-harness.log"
        with java_log.open("w") as log:
            subprocess.run(
                [
                    "java", f"-Xmx{args.java_heap}", "-Djava.awt.headless=true",
                    f"-Duser.home={java_home}",
                    "-cp", classpath_value,
                    "BaseTissueCropperReferenceHarness", str(raw_input), str(java_output),
                    str(shape[0]), str(shape[1]), str(shape[2]), str(crop_edge),
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )

        java_metadata = json.loads((java_output / "java-crop.json").read_text())
        java_bounds = tuple(tuple(int(v) for v in axis) for axis in java_metadata["bounds_inclusive"])

        from pymacruise.cropping import CropParameters, crop_image, find_crop_parameters

        python_crop_parameters = find_crop_parameters(
            gm_mask, threshold=0.0, edge_width=crop_edge,
            background=float(parameters.get("background", 0.0)),
        )
        python_bounds = python_crop_parameters.bounds
        if python_bounds != manifest_bounds:
            raise ValueError(
                f"Saved Atlas_GM crop bounds {manifest_bounds} do not match Python crop API {python_bounds}"
            )
        bounds_equal = java_bounds == python_bounds
        crop_shape = tuple(int(v) for v in java_metadata["cropped_shape"])

        crop_parameters = CropParameters(
            original_shape=shape,
            bounds=python_bounds,
            background=float(parameters.get("background", 0.0)),
            original_affine=python_affine,
        )
        expected_affine = python_affine.copy()
        translation = np.eye(4, dtype=np.float64)
        translation[:3, 3] = [axis[0] for axis in python_bounds]
        expected_affine = expected_affine @ translation

        results = []
        for name, (source_path, integer, saved_key) in targets.items():
            source = nifti(source_path)
            python_crop = crop_image(source, crop_parameters)
            if saved_key:
                saved_path = attempt_dir / f"{saved_key}.nii.gz"
                saved_crop = nifti(saved_path)
                assert_same_grid(python_crop, saved_crop, f"saved {saved_key}")
                python_data = np.asanyarray(saved_crop.dataobj)
                saved_path_record = record_input(saved_path)
            else:
                python_data = np.asanyarray(python_crop.dataobj)
                saved_path_record = None
            java_path = java_output / f"{name}{JAVA_SUFFIX[integer]}"
            java_data = java_image_data(java_path, integer)
            comparison = compare_arrays(java_data, python_data, name, crop_shape)
            comparison["java_output_bytes"] = java_path.stat().st_size
            comparison["python_crop_affine_matches_expected"] = bool(
                np.allclose(python_crop.affine, expected_affine, rtol=0.0, atol=1e-5)
            )
            comparison["saved_python_output"] = saved_path_record
            results.append(comparison)
        shutil.rmtree(java_output)

    report = {
        "comparison": "selected compiled JIST CubicVolumeCropper versus saved Python base-tissue-preparation crop",
        "run_root": str(run_root),
        "run_manifest": record_input(run_manifest),
        "attempt": attempt_number,
        "attempt_dir": str(attempt_dir),
        "attempt_provenance": record_input(provenance_path),
        "workflow_stage_record_status": stage_record.get("status"),
        "workflow_stage_record_attempt": stage_record.get("attempt"),
        "maseg_input": record_input(maseg_path),
        "bounding_mask": record_input(mask_path),
        "jist_runtime_jar": record_input(RUNTIME_JAR),
        "compiled_cropper_class_sha256": jar_entry_sha256(
            RUNTIME_JAR, "edu/jhmi/rad/medic/utilities/CubicVolumeCropper.class"
        ),
        "compiled_crop_parameters_class_sha256": jar_entry_sha256(
            RUNTIME_JAR, "edu/jhmi/rad/medic/utilities/CropParameters.class"
        ),
        "source_wrapper": str(WRAPPER_SOURCE),
        "source_wrapper_sha256": sha256(WRAPPER_SOURCE),
        "java_harness": str(HARNESS),
        "java_log": str(java_log),
        "java_heap_cap": args.java_heap,
        "source_wrapper": "CropBasedOnGM calls crop(gm, 0, edge), gets last params, then crop(target, params)",
        "threshold": 0.0,
        "edge_width": crop_edge,
        "input_shape": list(shape),
        "ras_axis_codes": list(nib.aff2axcodes(gm_mask.affine)),
        "python_bounds_inclusive": [list(axis) for axis in python_bounds],
        "java_bounds_inclusive": [list(axis) for axis in java_bounds],
        "bounds_match": bounds_equal,
        "cropped_shape": list(crop_shape),
        "python_crop_affine": expected_affine.tolist(),
        "scope": {
            "direct_class": "compiled JIST CubicVolumeCropper and CropParameters from the supplied selected runtime JAR",
            "compared_arrays": [item["name"] for item in results],
            "java_header_geometry": "not emitted by CubicVolumeCropper; Python output affine checked against offset formula separately",
            "source_separation": "inputs are saved Python-produced maps; Java cortical separation is not replayed here",
            "source_transfer": "uses already-transferred second-pass TOADS maps; no transfer is rerun here",
            "edge_limit": "this human MaSeg crop box is interior to the input image, so zero-padding edge behavior is not exercised",
            "historical_limit": "fixed-input crop-stage comparison only; not whole-workflow equivalence or anatomical acceptance",
        },
        "results": results,
        "all_arrays_equal": bool(
            bounds_equal
            and all(
                item.get("exact_values_match")
                and item.get("shape_match")
                and item.get("python_crop_affine_matches_expected")
                for item in results
            )
        ),
    }
    report_path = output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Java crop bounds: {java_bounds}")
    print(f"Python crop bounds: {python_bounds}")
    for result in results:
        print(f"{result['name']}: {result.get('different_voxels')} differing voxels")
    print(f"Report: {report_path}")
    if not report["all_arrays_equal"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
