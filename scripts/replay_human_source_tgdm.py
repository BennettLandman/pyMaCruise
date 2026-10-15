#!/usr/bin/env python3
"""Replay selected source GVF and TGDM on the downloaded human BASE inputs.

The exported human result contains only a scalar GVF NIfTI although the
selected source TGDM requires a three-component vector. This runner regenerates
the vector from the exported ACE-GM crop using the pinned VandyCRUISE Java
source, checks the saved scalar against the regenerated first component, then
executes the selected Java TGDM sequence against the exact exported inputs.
The result is a held-fixed source-stage replay, not a fresh full workflow.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("/Users/landmaba/dev/mipav/sample-macruise/output")
TOADS_SOURCE = Path("/Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26")
MATLAB_SOURCE = Path("/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE")
sys.path.insert(0, str(ROOT / "src"))

from pymacruise.cropping import find_crop_parameters  # noqa: E402
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics  # noqa: E402
from pymacruise.surfaces import extract_levelset_surface  # noqa: E402
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_arrays  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_record(relative: Path, original: Path | None = None) -> dict[str, object]:
    copied = ROOT / "tests/java_reference/upstream" / relative
    record: dict[str, object] = {"harness_source_copy": str(copied), "sha256": _sha256(copied)}
    if original is not None:
        original = original.resolve()
        original_hash = _sha256(original)
        record.update({
            "selected_source": str(original),
            "selected_source_sha256": original_hash,
            "copy_matches_selected_source": original_hash == record["sha256"],
        })
        if original_hash != record["sha256"]:
            raise ValueError(f"harness source copy differs from selected source: {relative}")
    return record


def _compile_gvf(classes: Path) -> None:
    source_root = ROOT / "tests/java_reference"
    upstream = source_root / "upstream"
    gvf_sources = sorted(
        path for path in (upstream / "edu/jhu/ece/iacl/algorithms/gvf").glob("*.java")
        if path.name != "FastMarchingGradient.java"
    )
    subprocess.run(
        [
            "javac", "-d", str(classes), "-sourcepath",
            f"{source_root / 'stubs'}:{upstream}",
            *(str(path) for path in gvf_sources),
            str(source_root / "GVFReferenceHarness.java"),
        ], check=True, cwd=ROOT,
    )


def _compile_tgdm(classes: Path) -> None:
    source_root = ROOT / "tests/java_reference"
    upstream = source_root / "upstream"
    subprocess.run(
        [
            "javac", "-d", str(classes), "-sourcepath",
            f"{source_root / 'stubs'}:{upstream}",
            str(source_root / "PlusTGDMReferenceHarness.java"),
        ], check=True, cwd=ROOT,
    )


def _write_scalar_input(path: Path, values: np.ndarray) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *values.shape))
        stream.write(np.asarray(values, dtype=">f4").tobytes(order="C"))


def _write_tgdm_input(path: Path, arrays: list[np.ndarray]) -> None:
    shape = arrays[0].shape
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for array in arrays:
            if array.shape[:3] != shape or array.ndim not in (3, 4):
                raise ValueError("all TGDM inputs must use one crop shape with optional vector channels")
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="C"))


def _read_gvf_output(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    values = path.read_bytes()
    dimensions = struct.unpack(">III", values[:12])
    if dimensions != shape:
        raise ValueError(f"Java GVF output shape {dimensions} != {shape}")
    vector = np.frombuffer(values, dtype=">f4", offset=12)
    if vector.size != int(np.prod(shape)) * 3:
        raise ValueError("Java GVF result does not contain three components per voxel")
    return vector.astype(np.float32).reshape((*shape, 3))


def _read_tgdm_output(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    with gzip.open(path, "rb") as stream:
        result = np.frombuffer(stream.read(), dtype=">f4")
    if result.size != int(np.prod(shape)):
        raise ValueError(f"Java TGDM output {path.name} has an unexpected size")
    return result.astype(np.float32).reshape(shape)


def _surface_metrics(java_values: np.ndarray, source_values: np.ndarray, affine: np.ndarray) -> dict:
    image_java = nib.Nifti1Image(java_values, affine)
    image_source = nib.Nifti1Image(source_values, affine)
    java_surface = extract_levelset_surface(image_java)
    source_surface = extract_levelset_surface(image_source)
    return {
        "dice_at_zero": float(dice_coefficient(java_values <= 0, source_values <= 0)),
        **mesh_surface_distance_metrics(
            java_surface.vertices, java_surface.faces,
            source_surface.vertices, source_surface.faces,
        ),
        "java_vertices": int(len(java_surface.vertices)),
        "source_vertices": int(len(source_surface.vertices)),
        "java_faces": int(len(java_surface.faces)),
        "source_faces": int(len(source_surface.faces)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-output", type=Path, default=SOURCE)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    source, work, report_path = args.source_output.resolve(), args.work_dir.resolve(), args.report.resolve()
    work.mkdir(parents=True, exist_ok=True)
    gvf_classes, tgdm_classes = work / "gvf-classes", work / "tgdm-classes"
    _compile_gvf(gvf_classes)
    _compile_tgdm(tgdm_classes)

    atlas_gm_path = source / "CorticalSeperation/Atlas_GM.nii.gz"
    atlas_gm = nib.load(str(atlas_gm_path))
    crop = find_crop_parameters(atlas_gm, threshold=0.0, edge_width=5)
    bounds = crop.bounds
    crop_slices = tuple(slice(lo, hi + 1) for lo, hi in bounds)
    shape = crop.cropped_shape

    tgdm_input_dir = source / "TGDMinputs"
    input_paths = {
        "membership": tgdm_input_dir / "target_image_Levelimg_uc_invcp_refspace.nii.gz",
        "white_matter": tgdm_input_dir / "target_image_WMimg_uc_invcp_refspace.nii.gz",
        "gray_matter": tgdm_input_dir / "target_image_GMimg_uc_invcp_refspace.nii.gz",
        "ace_gray_matter": tgdm_input_dir / "target_image_ACEGMimg_uc_invcp_refspace.nii.gz",
        "ace_skeleton": tgdm_input_dir / "target_image_ACEskelectonimg_uc_invcp_refspace.nii.gz",
        "saved_scalar_gvf": tgdm_input_dir / "target_image_GVFimg_uc_invcp_refspace.nii.gz",
    }
    arrays: dict[str, np.ndarray] = {}
    geometry = []
    for name, path in input_paths.items():
        image = nib.load(str(path))
        if image.ndim != 3 or image.shape != atlas_gm.shape:
            raise ValueError(f"{name} is not on the full source grid")
        values = np.asarray(image.dataobj, dtype=np.float32)
        if not np.isfinite(values).all():
            raise ValueError(f"{name} contains nonfinite values")
        arrays[name] = values[crop_slices].copy()
        geometry.append({"name": name, "shape": list(image.shape), "affine": image.affine.tolist(), "sha256": _sha256(path)})

    gvf_input = work / "gvf-crop.bin"
    gvf_output = work / "gvf-java.bin"
    _write_scalar_input(gvf_input, arrays["ace_gray_matter"])
    subprocess.run(
        ["java", "-Xms256m", "-Xmx4g", "-cp", str(gvf_classes),
         "GVFReferenceHarness", str(gvf_input), str(gvf_output), "1.0", "200"],
        check=True, cwd=ROOT,
    )
    vector = _read_gvf_output(gvf_output, shape)
    saved_gvf = arrays["saved_scalar_gvf"]
    vector_component_diagnostic = {}
    for component in range(3):
        difference = np.abs(vector[..., component].astype(np.float64) - saved_gvf.astype(np.float64))
        vector_component_diagnostic[str(component)] = {
            "mae_vs_exported_scalar": float(difference.mean()),
            "max_abs_vs_exported_scalar": float(difference.max(initial=0.0)),
            "correlation_vs_exported_scalar": float(np.corrcoef(vector[..., component].ravel(), saved_gvf.ravel())[0, 1]),
        }

    tgdm_order = ("membership", "white_matter", "gray_matter", "ace_gray_matter", "ace_skeleton")
    tgdm_input = work / "tgdm-input.bin"
    _write_tgdm_input(tgdm_input, [arrays[name] for name in tgdm_order] + [vector])
    java_output = work / "java-tgdm"
    java_output.mkdir(exist_ok=True)
    subprocess.run(
        ["java", "-Xms256m", "-Xmx6g", "-cp", str(tgdm_classes),
         "PlusTGDMReferenceHarness", str(tgdm_input), str(java_output)],
        check=True, cwd=ROOT,
    )

    # Run the Python array implementation on exactly the same crop and vector
    # as the Java harness. This isolates TGDM implementation parity from NIfTI
    # orientation, crop construction, GVF generation, and upstream workflows.
    python_execution = tgdm_nested_arrays(
        arrays["membership"], arrays["white_matter"], arrays["gray_matter"],
        arrays["ace_gray_matter"], ace_skeleton=arrays["ace_skeleton"],
        gvf=vector,
        parameters=TGDMBaselineParameters(white_matter_iso=0.5),
    )
    python_products = {
        "central": python_execution.central.level_set,
        "inner": python_execution.inner.level_set,
        "outer": python_execution.outer.level_set,
    }

    raw_t1 = nib.load(str(source / "Raw/T1.nii.gz"))
    cropped_affine = raw_t1.affine.copy()
    translation = np.eye(4)
    translation[:3, 3] = crop.offset
    cropped_affine = cropped_affine @ translation
    source_products = {
        "central": source / "Surfaces/target_image_GMimg_centralVol.nii.gz",
        "inner": source / "Surfaces/target_image_GMimg_innerVol.nii.gz",
        "outer": source / "Surfaces/target_image_GMimg_outerVol.nii.gz",
    }
    comparisons = {}
    for name, path in source_products.items():
        source_image = nib.load(str(path))
        if source_image.shape != atlas_gm.shape:
            raise ValueError(f"saved source product {name} has an unexpected grid")
        source_crop = np.asarray(source_image.dataobj, dtype=np.float32)[crop_slices]
        java_values = _read_tgdm_output(java_output / f"{name}.raw.gz", shape)
        python_values = np.asarray(python_products[name], dtype=np.float32)
        difference = np.abs(java_values.astype(np.float64) - source_crop.astype(np.float64))
        python_java_difference = np.abs(python_values.astype(np.float64) - java_values.astype(np.float64))
        comparisons[name] = {
            "differing_voxels": int(np.count_nonzero(java_values != source_crop)),
            "sign_differing_voxels": int(np.count_nonzero(np.signbit(java_values) != np.signbit(source_crop))),
            "mean_absolute_voxel_error": float(difference.mean()),
            "maximum_absolute_voxel_error": float(difference.max(initial=0.0)),
            "surface_metrics": _surface_metrics(java_values, source_crop, cropped_affine),
            "python_vs_java": {
                "differing_voxels": int(np.count_nonzero(python_values != java_values)),
                "sign_differing_voxels": int(np.count_nonzero(np.signbit(python_values) != np.signbit(java_values))),
                "mean_absolute_voxel_error": float(python_java_difference.mean()),
                "maximum_absolute_voxel_error": float(python_java_difference.max(initial=0.0)),
                "surface_metrics": _surface_metrics(python_values, java_values, cropped_affine),
            },
        }

    result = {
        "case": "downloaded human BASE source TGDM replay",
        "evidence_category": "held-fixed Java source-stage replay using exported historical workflow inputs",
        "scope_limit": "not a fresh full MaCRUISE Java/MATLAB run; the exported scalar GVF is checked as a missing-vector diagnostic and the three-component vector is regenerated from the exported ACE-GM crop",
        "coordinate_system": "RAS millimetres for surface distances",
        "source_output": str(source),
        "source_atlas_gm_sha256": _sha256(atlas_gm_path),
        "source_output_t1_sha256": _sha256(source / "Raw/T1.nii.gz"),
        "source_output_segmentation_sha256": _sha256(source / "Raw/Seg.nii.gz"),
        "full_shape": list(atlas_gm.shape),
        "crop_bounds_inclusive": [list(pair) for pair in bounds],
        "crop_shape": list(shape),
        "crop_affine_ras": cropped_affine.tolist(),
        "java_input_hashes_and_geometry": geometry,
        "saved_scalar_gvf_component_diagnostic": vector_component_diagnostic,
        "source_surface_comparisons": comparisons,
        "java_source": {
            "gvf": str(ROOT / "tests/java_reference/GVFReferenceHarness.java"),
            "tgdm": str(ROOT / "tests/java_reference/PlusTGDMReferenceHarness.java"),
        },
        "source_file_hashes": {
            "ma_generic_tgdm": _source_record(
                Path("edu/vanderbilt/masi/plugins/CRUISE/tgdm/MaGenericTGDM.java"),
                TOADS_SOURCE / "toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/tgdm/MaGenericTGDM.java",
            ),
            "distance_field": _source_record(
                Path("edu/vanderbilt/masi/algorithms/CRUISE/ace/DistanceFieldyh.java"),
                TOADS_SOURCE / "toads-cruise/src/edu/vanderbilt/masi/algorithms/CRUISE/ace/DistanceFieldyh.java",
            ),
            "generalized_gvf": _source_record(
                Path("edu/jhu/ece/iacl/algorithms/gvf/GradVecFlowGeneralized.java"),
                TOADS_SOURCE / "VandyCRUISE/src/edu/jhu/ece/iacl/algorithms/gvf/GradVecFlowGeneralized.java",
            ),
            "ma_tgdm_crop_wrapper": {
                "path": str(TOADS_SOURCE / "toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/tgdm/MaTGDMcrop.java"),
                "sha256": _sha256(TOADS_SOURCE / "toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/tgdm/MaTGDMcrop.java"),
                "used_as": "parameter and stage-order source; the Java CLI wrapper itself is not run",
            },
            "matlab_cli_wrapper": {
                "path": str(MATLAB_SOURCE / "pipeline/run_TGDMinvcrop_jist.m"),
                "sha256": _sha256(MATLAB_SOURCE / "pipeline/run_TGDMinvcrop_jist.m"),
                "used_as": "selected workflow parameter defaults and command input ordering",
            },
        },
        "source_selected_parameters": {
            "gvf_lambda": 1.0,
            "gvf_max_iterations": 200,
            "connectivity": "18/6",
            "initial_iso": 126.9,
            "wm_iso": 0.5,
            "gm_iso": 0.4,
            "initial_smoothing_iterations": 6,
            "inner_iterations": 4,
            "central_iterations": 7,
            "outer_iterations": 7,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
