#!/usr/bin/env python3
"""Compare the committed fine-grid LesionTOADS Dura fixture with selected Java.

The comparison uses only committed `.nii.gz` inputs and a pinned Java excerpt;
it does not need the ignored workflow build or the original source checkout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile

import nibabel as nib
import numpy as np


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests/data/synthetic/plus_fine_grid_dura_iteration_6"
HARNESS = REPO / "tests/java_reference/LesionToadsDuraCheckpointReferenceHarness.java"
PINNED_JAVA_SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def checked_path(record: dict, fixture: Path, label: str) -> Path:
    path = (fixture / record["path"]).resolve()
    if not path.is_file() or path.is_symlink() or digest(path) != record.get("sha256"):
        raise ValueError(f"{label} missing or failed saved SHA-256 check: {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-dir", type=Path, default=FIXTURE)
    parser.add_argument("--iteration", type=int, default=6,
                        help="expected completed iteration from the saved loop summary")
    parser.add_argument("--report", type=Path,
                        help="machine-readable result (default: build/p41z-dura-java-comparison.json)")
    args = parser.parse_args()
    fixture = args.fixture_dir.resolve()
    provenance_path = fixture / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if not provenance.get("synthetic_only") or provenance.get("spatial_convention") != "RAS":
        raise ValueError("fixture provenance must mark synthetic RAS data")
    completed_iteration = int(provenance["completed_iteration"])
    if completed_iteration != args.iteration:
        raise ValueError(f"expected completed iteration {args.iteration}, got {completed_iteration}")
    class_names = provenance["class_names"]
    if len(class_names) != 13:
        raise ValueError(f"expected 13 atlas classes, got {len(class_names)}")
    memberships_path = checked_path(provenance["files"]["memberships_15ch"], fixture, "saved 15-channel memberships")
    segmentation_path = checked_path(provenance["files"]["segmentation"], fixture, "saved segmentation")
    memberships_img = nib.load(str(memberships_path))
    segmentation_img = nib.load(str(segmentation_path))
    if memberships_img.shape != (*segmentation_img.shape, 15):
        raise ValueError("saved membership fixture is not xyz-by-15")
    if not np.allclose(memberships_img.affine, segmentation_img.affine, rtol=0, atol=1e-5):
        raise ValueError("saved membership and segmentation grids differ")
    memberships = np.asarray(memberships_img.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_img.dataobj, dtype=np.int32)
    lesion = memberships[..., 13]
    black_hole = memberships[..., 14]
    classes = memberships[..., :13]
    shape = segmentation.shape
    if memberships.shape[:3] != shape:
        raise ValueError("saved arrays have inconsistent dimensions")

    output_paths = {
        name: checked_path(provenance["files"][f"python_{name}"], fixture, f"saved Python {name}")
        for name in ("csf", "gm", "wm_fill", "wm_mask")
    }
    python_outputs = []
    for name in ("csf", "gm", "wm_fill", "wm_mask"):
        image = nib.load(str(output_paths[name]))
        if image.shape != shape or not np.allclose(image.affine, segmentation_img.affine, rtol=0, atol=1e-5):
            raise ValueError(f"saved Python {name} output is on a different grid")
        python_outputs.append(np.asarray(image.dataobj, dtype=np.float32))

    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("JDK javac and java are required")
    if provenance.get("selected_java_source_sha256") != PINNED_JAVA_SOURCE_SHA256:
        raise ValueError("fixture does not identify the pinned selected Java source method")

    report_path = args.report or REPO / "build/p41z-dura-java-comparison.json"
    report_path = report_path.resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pymacruise-p41z-dura-") as temporary:
        temp = Path(temporary)
        classes_dir = temp / "classes"
        classes_dir.mkdir()
        subprocess.run([javac, "-d", str(classes_dir), str(HARNESS)], check=True, capture_output=True, text=True)
        input_path, reference_path = temp / "input.bin", temp / "java-outputs.bin"
        with input_path.open("wb") as out:
            out.write(struct.pack(">IIII", *shape, len(class_names)))
            out.write(classes.astype(">f4", copy=False).tobytes(order="C"))
            out.write(segmentation.astype(">i4", copy=False).tobytes(order="C"))
            out.write(lesion.astype(">f4", copy=False).tobytes(order="C"))
            out.write(black_hole.astype(">f4", copy=False).tobytes(order="C"))
        subprocess.run(
            [java, "-cp", str(classes_dir), HARNESS.stem, str(input_path), str(reference_path), *class_names],
            check=True, capture_output=True, text=True,
        )
        java_output = np.fromfile(reference_path, dtype=">f4").astype(np.float32)
    expected_size = 4 * int(np.prod(shape))
    if java_output.size != expected_size:
        raise ValueError(f"Java returned {java_output.size} values, expected {expected_size}")
    java_outputs = java_output.reshape((4, *shape))
    channel_names = ("csf", "gm", "wm_fill", "wm_mask")
    comparisons = {}
    for index, name in enumerate(channel_names):
        actual = java_outputs[index]
        reference = python_outputs[index]
        mismatch = int(np.count_nonzero(actual != reference))
        comparisons[name] = {
            "exact_float32_match": mismatch == 0,
            "different_voxels": mismatch,
            "voxel_count": int(actual.size),
            "max_abs_difference": float(np.max(np.abs(actual - reference), initial=0.0)),
        }
    result = {
        "case": provenance["case"],
        "completed_iteration": completed_iteration,
        "synthetic_only": True,
        "shape_xyz": list(shape),
        "class_names": class_names,
        "memberships_used": "channels 0-12",
        "lesion_membership_channel": 13,
        "black_hole_membership_channel": 14,
        "fixture_provenance_sha256": digest(provenance_path),
        "source_checkpoint": provenance["source_checkpoint"],
        "java_source": provenance["selected_java_source"],
        "java_source_sha256": provenance["selected_java_source_sha256"],
        "java_harness_sha256": digest(HARNESS),
        "inputs": {
            "memberships": digest(memberships_path),
            "segmentation": digest(segmentation_path),
        },
        "python_output_hashes": {name: digest(path) for name, path in output_paths.items()},
        "comparisons": comparisons,
        "all_four_exact": all(item["exact_float32_match"] for item in comparisons.values()),
    }
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_four_exact"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
