#!/usr/bin/env python3
"""Compare the source-selected JIST transfer with Python on saved human inputs.

This uses the selected VandyCRUISE ReorientVolume source with the existing
compiled MIPAV/JIST runtime artifacts. Java output is compared voxel-for-voxel
to ``pymacruise.spatial_transfer.transfer_to_reference_space``. The three
source maps are processed sequentially to keep peak memory bounded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np


REPO = Path(__file__).resolve().parents[1]
MIPAV = REPO.parent / "mipav-src"
JIST = REPO.parent / "jist-scm-2026-09-26"
VANDY = REPO.parent / "toads-cruise-scm-2026-09-26" / "toads-cruise"
BUILD = REPO / "build" / "java-reference" / "reorient-volume-human-transfer"
HARNESS = REPO / "tests" / "java_reference" / "ReorientVolumeHumanTransferReferenceHarness.java"
REORIENT = VANDY / "src" / "edu" / "jhu" / "ece" / "iacl" / "algorithms" / "volume" / "ReorientVolume.java"

DEFAULT_RUN = REPO / "build" / "invivo-base-2026-10-10-ras-grid"
DEFAULT_INPUT_ROOT = REPO / "build" / "invivo-base-2026-10-10-inputs" / "ras-grid"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classpath() -> str:
    entries = [
        str(BUILD / "classes"),
        str(REPO / "build" / "java-reference" / "vandycruise-core-classes"),
        str(REPO / "build" / "java-reference" / "toads-orientation-mapping-classes"),
        str(REPO / "build" / "java-reference" / "clean-source-compile"),
        str(REPO / "build" / "java-reference" / "jist-runtime-without-selected-conflicts.jar"),
    ]
    # Include only supplied runtime libraries. Compiled MIPAV classes are first
    # so similarly named classes in legacy application jars cannot override them.
    for library_root in (MIPAV / "src" / "lib", MIPAV / "src" / "libnative", JIST / "JIST" / "lib"):
        entries.extend(str(path) for path in sorted(library_root.rglob("*.jar")))
    return os.pathsep.join(entries)


def ensure_compiled() -> str:
    needed = (REORIENT, HARNESS)
    missing = [str(path) for path in needed if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing selected Java source(s): " + ", ".join(missing))
    for required in (
        REPO / "build/java-reference/vandycruise-core-classes/gov/nih/mipav/model/algorithms/AlgorithmTransform.class",
        REPO / "build/java-reference/vandycruise-core-classes/gov/nih/mipav/model/file/FileIO.class",
        REPO / "build/java-reference/toads-orientation-mapping-classes/edu/jhmi/rad/medic/structures/RotationMatrix.class",
        REPO / "build/java-reference/jist-runtime-without-selected-conflicts.jar",
    ):
        if not required.exists():
            raise FileNotFoundError(f"Required precompiled Java runtime artifact is missing: {required}")
    classes = BUILD / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    generated_source = BUILD / "generated-src" / "edu" / "jhu" / "ece" / "iacl" / "algorithms" / "volume" / "ReorientVolume.java"
    generated_source.parent.mkdir(parents=True, exist_ok=True)
    source_text = REORIENT.read_text()
    # The selected 2014 source calls the old uppercase TransMatrix API. The
    # supplied precompiled MIPAV runtime exposes the same operations as
    # WildMagic Matrix4f.identity()/set(). Rewrite only these matrix API names
    # in a build-local copy; the orientation, dimensions, interpolation, and
    # AlgorithmTransform call remain the selected source statements.
    if source_text.count("transform.MakeIdentity();") != 2 or source_text.count("transform.Set(") != 4:
        raise RuntimeError("Selected ReorientVolume matrix API occurrences changed; refusing an implicit rewrite")
    generated_source.write_text(
        source_text.replace("transform.MakeIdentity();", "transform.identity();")
        .replace("transform.Set(", "transform.set(")
    )
    cp = classpath()
    subprocess.run(
        ["javac", "-classpath", cp, "-d", str(classes), str(generated_source), str(HARNESS)],
        check=True,
    )
    return cp


def image_geometry(image: nib.spatialimages.SpatialImage) -> dict[str, object]:
    return {
        "shape": list(image.shape),
        "axis_codes": list(nib.aff2axcodes(image.affine)),
        "affine": np.asarray(image.affine, dtype=np.float64).tolist(),
        "dtype": str(image.get_data_dtype()),
    }


def select_completed_toads_attempt(run_root: Path, attempt: int | None) -> tuple[int, Path]:
    """Select a finished TOADS attempt by its completion marker and outputs.

    A directory may exist while a retry is still running, so directory order
    alone is not enough. The stage log's final completion message is written
    only after TOADS has emitted its dura-input images.
    """
    stage_root = run_root / "stages" / "toads-second"
    if attempt is None:
        candidates: list[tuple[int, Path]] = []
        for path in stage_root.glob("attempt-*"):
            if not path.is_dir():
                continue
            suffix = path.name.removeprefix("attempt-")
            if suffix.isdigit():
                candidates.append((int(suffix), path))
        candidates.sort(reverse=True)
    else:
        candidates = [(attempt, stage_root / f"attempt-{attempt:02d}")]

    required = (
        "dura-inputs/stripped_t1_gm.nii.gz",
        "dura-inputs/stripped_t1_wmfill.nii.gz",
        "dura-inputs/stripped_t1_wmmask.nii.gz",
    )
    rejected: list[str] = []
    for number, path in candidates:
        log_path = path / "stage.log"
        if not log_path.is_file():
            rejected.append(f"{path}: missing stage.log")
            continue
        log_text = log_path.read_text(errors="replace")
        if "TOADS segmentation loop complete" not in log_text:
            rejected.append(f"{path}: no TOADS completion marker in stage.log")
            continue
        missing = [name for name in required if not (path / name).is_file()]
        if missing:
            rejected.append(f"{path}: missing required outputs: {', '.join(missing)}")
            continue
        return number, path

    if attempt is not None:
        reason = "; ".join(rejected) if rejected else "attempt directory does not exist"
        raise FileNotFoundError(f"TOADS attempt {attempt} is not complete under {stage_root}: {reason}")
    detail = "; ".join(rejected) if rejected else "no attempt-NN directories found"
    raise FileNotFoundError(f"No completed TOADS attempt found under {stage_root}: {detail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN)
    parser.add_argument(
        "--attempt", type=int,
        help="TOADS attempt number; by default, use the newest attempt with a completion marker and required outputs",
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument(
        "--output-dir", type=Path,
        help="report directory (default: build/java-reference/reorient-volume-human-transfer/human-attempt-NN)",
    )
    args = parser.parse_args()

    selected_attempt, attempt_root = select_completed_toads_attempt(args.run_root, args.attempt)
    toads_root = attempt_root / "dura-inputs"
    if args.output_dir is None:
        args.output_dir = BUILD / f"human-attempt-{selected_attempt:02d}"
    reference = args.input_root / "orig_target_seg_ras_uint16.nii.gz"
    sources = {
        "gm": toads_root / "stripped_t1_gm.nii.gz",
        "wmfill": toads_root / "stripped_t1_wmfill.nii.gz",
        "wmmask": toads_root / "stripped_t1_wmmask.nii.gz",
    }
    for path in (reference, *sources.values()):
        if not path.is_file():
            raise FileNotFoundError(f"Human transfer input is missing: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    java_home = BUILD / "java-home"
    java_home.mkdir(parents=True, exist_ok=True)

    cp = ensure_compiled()
    java_args = [
        "java", "-Djava.awt.headless=true", f"-Duser.home={java_home}",
        "-Xms128m", "-Xmx2g", "-cp", cp,
        "ReorientVolumeHumanTransferReferenceHarness", str(reference),
    ]
    for label, source in sources.items():
        stem = args.output_dir / label
        java_args.extend((label, str(source), f"{stem}.java.f32be,{stem}.java.json"))
    java_log = args.output_dir / "java-harness.log"
    with java_log.open("w") as log:
        completed = subprocess.run(
            java_args, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    if completed.returncode:
        tail = "\n".join(java_log.read_text(errors="replace").splitlines()[-35:])
        raise RuntimeError(f"Java transfer harness exited {completed.returncode}; tail of {java_log}:\n{tail}")

    # Import late so compiling/running Java works even if package import setup is
    # supplied externally through PYTHONPATH.
    from pymacruise.spatial_transfer import transfer_to_reference_space

    reference_img = nib.load(str(reference))
    report: dict[str, object] = {
        "comparison": "selected JIST ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR) versus Python spatial_transfer.transfer_to_reference_space",
        "run_root": str(args.run_root),
        "selected_toads_attempt": selected_attempt,
        "toads_attempt_root": str(attempt_root),
        "reference": {"path": str(reference), "sha256": sha256(reference), "geometry": image_geometry(reference_img)},
        "java_source": str(REORIENT),
        "java_source_sha256": sha256(REORIENT),
        "java_harness": str(HARNESS),
        "mipav_runtime_classes": str(REPO / "build/java-reference/vandycruise-core-classes"),
        "jist_runtime_jar": str(REPO / "build/java-reference/jist-runtime-without-selected-conflicts.jar"),
        "compile_compatibility_rewrites": [
            "build-local ReorientVolume copy maps two legacy TransMatrix.MakeIdentity() calls to the supplied WildMagic Matrix4f.identity() API",
            "build-local ReorientVolume copy maps four legacy TransMatrix.Set(row,col,value) calls to the supplied Matrix4f.set(row,col,value) API",
        ],
        "results": {},
    }
    results = report["results"]
    assert isinstance(results, dict)
    for label, source in sources.items():
        source_img = nib.load(str(source))
        python_img = transfer_to_reference_space(source_img, reference_img, cval=0.0)
        expected = np.asarray(python_img.dataobj).astype(np.float32, copy=False)
        raw_path = args.output_dir / f"{label}.java.f32be"
        actual_linear = np.fromfile(raw_path, dtype=">f4")
        java_geometry_path = args.output_dir / f"{label}.java.json"
        java_geometry = json.loads(java_geometry_path.read_text())
        java_axis_orientation = java_geometry["output_model_geometry_after_reference_header"]["axis_orientation"]
        # MIPAV stores this RAS input with x=L2R, y=A2P, z=I2S. Its in-memory
        # y index therefore runs opposite to NIfTI's +A index; undo that single
        # storage-axis reversal before comparing physical voxel coordinates.
        if nib.aff2axcodes(source_img.affine) != ("R", "A", "S") or java_axis_orientation != [2, 4, 5]:
            raise ValueError(
                f"Unsupported MIPAV/NIfTI axis mapping for {label}: "
                f"RAS={nib.aff2axcodes(source_img.affine)}, MIPAV={java_axis_orientation}"
            )
        expected_in_mipav_order = np.flip(expected, axis=1)
        expected_linear = expected_in_mipav_order.ravel(order="F")
        if actual_linear.size != expected_linear.size:
            raise ValueError(
                f"{label}: Java returned {actual_linear.size} values; Python grid has {expected_linear.size}"
            )
        actual_mipav = actual_linear.astype(np.float32, copy=False).reshape(expected.shape, order="F")
        actual = np.flip(actual_mipav, axis=1)
        delta = np.abs(actual - expected)
        results[label] = {
            "source": {"path": str(source), "sha256": sha256(source), "geometry": image_geometry(source_img)},
            "java_metadata": java_geometry,
            "python_output_geometry": image_geometry(python_img),
            "java_modelimage_to_nifti_ras_axis_mapping": "flip y axis (MIPAV A2P storage versus NIfTI +A array order)",
            "java_voxel_count": int(actual.size),
            "bitwise_equal_float32": bool(np.array_equal(actual, expected)),
            "different_voxels": int(np.count_nonzero(actual != expected)),
            "max_absolute_error": float(np.max(delta, initial=0.0)),
            "mean_absolute_error": float(np.mean(delta, dtype=np.float64)),
            "java_output_sha256": sha256(raw_path),
        }
        del source_img, python_img, expected, expected_linear, actual_linear, actual, delta

    report["all_three_bitwise_equal"] = all(
        item["bitwise_equal_float32"] for item in results.values()
    )
    report["scope_note"] = (
        "This is an actual-human-input comparison. The saved human source and reference grids are both RAS and aligned; it does not establish transfer behavior for changed orientation, origin, spacing, or shape."
    )
    report_path = args.output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["all_three_bitwise_equal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
