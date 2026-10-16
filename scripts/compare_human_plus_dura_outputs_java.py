#!/usr/bin/env python3
"""Compare saved human PLUS Dura maps with the selected Java output method."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests/java_reference/LesionToadsHumanDuraOutputsReferenceHarness.java"
SOURCE = Path(
    "/Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/"
    "VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
)
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
)
OUTPUT_NAMES = ("csf", "gm", "wm_fill", "wm_mask")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workflow-root", type=Path,
        default=ROOT / "build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    workflow = args.workflow_root.resolve()
    work = args.output_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    registered = workflow / "lesion_workflow/registered-loop"
    memberships_path = registered / "registered_memberships_17ch.nii.gz"
    labels_path = registered / "registered_segmentation.nii.gz"
    memberships_image = nib.load(str(memberships_path))
    labels_image = nib.load(str(labels_path))
    memberships = np.asarray(memberships_image.dataobj, dtype=np.float32)
    labels = np.asarray(labels_image.dataobj)
    shape = memberships_image.shape[:3]
    if memberships_image.shape != (*shape, 17) or labels.shape != shape:
        raise ValueError("registered PLUS arrays must be [x,y,z,17] memberships and [x,y,z] labels")
    if not np.array_equal(memberships_image.affine, labels_image.affine):
        raise ValueError("registered membership and segmentation affines differ")
    if tuple(str(code) for code in nib.aff2axcodes(memberships_image.affine)) != ("R", "A", "S"):
        raise ValueError("registered PLUS inputs must be RAS")
    if not np.isfinite(memberships).all() or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("membership maps must be finite and segmentation labels integer-valued")
    if int(labels.min()) < 0 or int(labels.max()) >= len(CLASS_NAMES):
        raise ValueError("segmentation contains a class index outside the selected 13-class atlas")

    java_input = work / "human-plus-dura-input.bin"
    with java_input.open("wb") as stream:
        stream.write(struct.pack(">4i", *shape, 13))
        for class_index in range(13):
            stream.write(np.asarray(memberships[..., class_index], dtype=">f4").tobytes(order="C"))
        stream.write(np.asarray(labels, dtype=np.uint8).tobytes(order="C"))
        stream.write(np.asarray(memberships[..., 13], dtype=">f4").tobytes(order="C"))
        stream.write(np.asarray(memberships[..., 14], dtype=">f4").tobytes(order="C"))

    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("a JDK is required for this selected-source comparison")
    classes = work / "classes"
    classes.mkdir(exist_ok=True)
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, cwd=ROOT)
    java_output = work / "java-dura-output.bin"
    subprocess.run(
        [java, "-Xms256m", "-Xmx2g", "-cp", str(classes),
         "LesionToadsHumanDuraOutputsReferenceHarness", str(java_input), str(java_output)],
        check=True, cwd=ROOT,
    )
    n_voxels = int(np.prod(shape))
    java_products = np.fromfile(java_output, dtype=">f4").astype(np.float32)
    if java_products.size != 4 * n_voxels:
        raise ValueError("Java output has an unexpected element count")
    java_products = java_products.reshape((4, *shape))

    expected_paths = {
        name: workflow / f"dura-outputs/lesion_toads_{name}.nii.gz"
        for name in OUTPUT_NAMES
    }
    comparisons = {}
    for channel, name in enumerate(OUTPUT_NAMES):
        expected_image = nib.load(str(expected_paths[name]))
        if expected_image.shape != shape or not np.array_equal(expected_image.affine, memberships_image.affine):
            raise ValueError(f"saved {name} output is not on the exact registered RAS grid")
        expected = np.asarray(expected_image.dataobj, dtype=np.float32)
        scaled = np.multiply(java_products[channel], np.float32(255.0), dtype=np.float32)
        difference = np.abs(scaled.astype(np.float64) - expected.astype(np.float64))
        comparisons[name] = {
            "differing_values_after_wrapper_scale": int(np.count_nonzero(scaled != expected)),
            "mean_absolute_scaled_error": float(difference.mean()),
            "maximum_absolute_scaled_error": float(difference.max(initial=0.0)),
            "saved_nonzero_voxels": int(np.count_nonzero(expected)),
            "java_nonzero_voxels": int(np.count_nonzero(scaled)),
        }

    source_sha = sha256(SOURCE) if SOURCE.is_file() else None
    if source_sha is not None and source_sha != SOURCE_SHA256:
        raise ValueError("local selected LesionToadSegmentation source no longer matches the pinned source hash")
    result = {
        "case": "human PLUS case 100260, registered LesionTOADS-to-Dura handoff",
        "evidence_category": "held-fixed selected Java method replay on Python-generated human memberships",
        "scope_limit": "does not compare the full LesionTOADS iteration or independent original PLUS workflow",
        "workflow_root": str(workflow),
        "shape_xyz": list(shape),
        "spatial_convention": "RAS",
        "atlas_class_names": list(CLASS_NAMES),
        "inputs": {
            "registered_memberships_sha256": sha256(memberships_path),
            "registered_segmentation_sha256": sha256(labels_path),
            "packed_java_input_sha256": sha256(java_input),
        },
        "outputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in expected_paths.items()
        },
        "selected_source": {
            "path": str(SOURCE), "sha256": source_sha or SOURCE_SHA256,
            "expected_sha256": SOURCE_SHA256,
            "harness": str(HARNESS), "harness_sha256": sha256(HARNESS),
        },
        "wrapper_scale": 255.0,
        "comparison": comparisons,
        "all_four_outputs_bit_exact_after_wrapper_scale": all(
            item["differing_values_after_wrapper_scale"] == 0 for item in comparisons.values()
        ),
    }
    report = args.report.resolve()
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
