#!/usr/bin/env python3
"""Compare saved base preparation crops with the pinned Java identity formula.

This uses ``PlusLesionTransferCropRasIdentityReferenceHarness``. It is a
formula reference for unchanged RAS grids, not an execution of JIST's
ReorientVolume or CubicVolumeCropper classes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


HARNESS_RELATIVE = Path(
    "tests/java_reference/PlusLesionTransferCropRasIdentityReferenceHarness.java"
)
HARNESS_SHA256 = "a3fb7afa0365609e513d0cb480e169643637144e3010909bd0409befec266c94"
EDGE_WIDTH = 5
FULL_KEYS = (
    "atlas_gm",
    "atlas_wm",
    "atlas_mask",
    "atlas_brain_tissue_mask",
    "toads_gm_transferred",
    "toads_wmfill_transferred",
    "toads_mask_transferred",
)
CROPPED_KEYS = (
    "atlas_gm_cropped",
    "atlas_wm_cropped",
    "atlas_mask_cropped",
    "atlas_brain_tissue_mask_cropped",
    "toads_gm_cropped",
    "toads_wmfill_cropped",
    "toads_mask_cropped",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_ras(path: str | Path, name: str) -> nib.spatialimages.SpatialImage:
    image = nib.load(str(path))
    return normalize_image_to_ras(image, source_name=name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    prep = args.preparation_dir.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    provenance_path = prep / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    parameters = provenance["parameters"]
    if parameters.get("remove_cerebellum") is not False:
        raise SystemExit("Expected corrected V3 preparation with remove_cerebellum=false")
    crop = parameters["crop"]
    bounds = tuple(tuple(int(v) for v in pair) for pair in crop["bounds_inclusive"])
    if crop.get("edge_width_voxels") != EDGE_WIDTH:
        raise SystemExit(f"Expected source-selected crop edge width {EDGE_WIDTH}")

    full = [
        _load_ras(provenance["outputs"][key], key)
        for key in FULL_KEYS
    ]
    labels_path = Path(provenance["inputs"]["maseg"])
    labels = _load_ras(labels_path, "MaSeg labels")
    all_images = [*full, labels]
    shape = tuple(int(size) for size in full[0].shape)
    if len(shape) != 3 or any(tuple(image.shape) != shape for image in all_images):
        raise SystemExit("Full-grid inputs do not share one 3D shape")
    affine = full[0].affine
    if any(not np.allclose(image.affine, affine, rtol=0, atol=1e-5) for image in all_images):
        raise SystemExit("Full-grid inputs do not share one physical RAS grid")
    if tuple(shape) != tuple(crop["original_shape"]):
        raise SystemExit("Full-grid image shape differs from crop provenance")

    slices = tuple(slice(low, high + 1) for low, high in bounds)
    # Confirm each saved TOADS crop is an exact voxel copy from its already-
    # transferred full-grid input. This does not replay the earlier resampling.
    same_grid_rows: dict[str, dict[str, object]] = {}
    for full_key, crop_key, image in zip(
        FULL_KEYS[4:], CROPPED_KEYS[4:], full[4:], strict=True
    ):
        source_crop = np.asarray(image.dataobj)[slices]
        saved_crop = np.asarray(
            _load_ras(provenance["outputs"][crop_key], crop_key).dataobj
        )
        same_grid_rows[full_key] = {
            "cropped_output": crop_key,
            "exact_voxel_copy_from_same_grid_input": bool(
                np.array_equal(source_crop, saved_crop)
            ),
            "voxel_count": int(source_crop.size),
            "changed_voxels": int(np.count_nonzero(source_crop != saved_crop)),
        }

    arrays = [np.asarray(image.dataobj, dtype=np.float32) for image in all_images]
    input_path = out / "crop-inputs.bin"
    with input_path.open("wb") as stream:
        stream.write(struct.pack(">4i", *shape, EDGE_WIDTH))
        for array in arrays:
            np.asarray(array, dtype=">f4").ravel(order="C").tofile(stream)

    harness = repo / HARNESS_RELATIVE
    harness_sha = _sha256(harness)
    if harness_sha != HARNESS_SHA256:
        raise SystemExit("Pinned Java formula harness source hash changed")
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise SystemExit("javac and java are required for the pinned formula reference")
    classes = out / "classes"
    classes.mkdir(exist_ok=True)
    subprocess.run([javac, "-d", str(classes), str(harness)], check=True)
    watchdog = repo / "scripts/run_with_rss_watchdog.py"
    output_path = out / "crop-reference.bin"
    java_result = subprocess.run(
        [
            sys.executable,
            str(watchdog),
            "--limit-gb",
            "2",
            "--interval-seconds",
            "0.2",
            "--report-interval-seconds",
            "15",
            "--",
            java,
            "-Xms64m",
            "-Xmx768m",
            "-cp",
            str(classes),
            "PlusLesionTransferCropRasIdentityReferenceHarness",
            str(input_path),
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with output_path.open("rb") as stream:
        low_high = struct.unpack(">6i", stream.read(24))
        java_shape = struct.unpack(">3i", stream.read(12))
        count = int(np.prod(java_shape))
        java_arrays = [
            np.frombuffer(stream.read(count * 4), dtype=">f4")
            .astype(np.float32)
            .reshape(java_shape)
            for _ in arrays
        ]
        if stream.read():
            raise SystemExit("Unexpected trailing bytes in Java formula output")

    java_bounds = tuple(
        (int(low_high[2 * axis]), int(low_high[2 * axis + 1]))
        for axis in range(3)
    )
    if java_bounds != bounds:
        raise SystemExit(f"Java formula bounds {java_bounds} differ from saved {bounds}")
    if tuple(java_shape) != tuple(high - low + 1 for low, high in bounds):
        raise SystemExit("Java formula output shape differs from inclusive crop bounds")

    rows: dict[str, dict[str, object]] = {}
    crop_image_paths = [provenance["outputs"][key] for key in CROPPED_KEYS]
    for index, key in enumerate(CROPPED_KEYS):
        saved = _load_ras(crop_image_paths[index], key)
        java_array = java_arrays[index]
        python_array = np.asarray(saved.dataobj, dtype=np.float32)
        if tuple(saved.shape) != tuple(java_shape):
            raise SystemExit(f"Saved crop {key} has unexpected shape {saved.shape}")
        rows[key] = {
            "exact": bool(np.array_equal(java_array, python_array)),
            "changed_voxels": int(np.count_nonzero(java_array != python_array)),
            "voxel_count": int(java_array.size),
            "maximum_absolute_error": float(
                np.abs(java_array.astype(np.float64) - python_array.astype(np.float64)).max(initial=0)
            ),
        }

    # MaSeg labels are the eighth held-fixed volume. The base preparation
    # exporter does not persist a cropped label image, so compare the Java
    # result against the corresponding direct full-grid crop slices.
    label_expected = arrays[7][slices]
    label_java = java_arrays[7]
    rows["maseg_labels_crop"] = {
        "exact": bool(np.array_equal(label_java, label_expected)),
        "changed_voxels": int(np.count_nonzero(label_java != label_expected)),
        "voxel_count": int(label_java.size),
        "maximum_absolute_error": float(
            np.abs(label_java.astype(np.float64) - label_expected.astype(np.float64)).max(initial=0)
        ),
        "reference": "direct MaSeg full-grid slice using saved inclusive bounds",
    }

    report = {
        "comparison_type": "held-fixed RAS identity-grid Java formula reference; not JIST class execution",
        "harness": str(harness),
        "harness_sha256": harness_sha,
        "preparation_provenance": str(provenance_path),
        "preparation_provenance_sha256": _sha256(provenance_path),
        "remove_cerebellum": False,
        "edge_width_voxels": EDGE_WIDTH,
        "full_grid_shape_xyz": list(shape),
        "ras_affine": affine.tolist(),
        "bounds_inclusive": [list(pair) for pair in bounds],
        "crop_shape_xyz": list(java_shape),
        "inputs": {
            key: {
                "path": provenance["outputs"][key],
                "sha256": _sha256(Path(provenance["outputs"][key])),
            }
            for key in FULL_KEYS
        },
        "labels_input": {"path": str(labels_path), "sha256": _sha256(labels_path)},
        "packed_input_sha256": _sha256(input_path),
        "same_grid_toads_identity_checks": same_grid_rows,
        "crop_outputs": rows,
        "all_crop_outputs_exact": all(bool(row["exact"]) for row in rows.values()),
        "watchdog_output": java_result.stdout.strip(),
        "scope_limitations": [
            "The pinned Java harness reproduces the identity-grid bounds, zero-padding, and voxel-copy formula; it does not invoke JIST CropBasedOnGM or CubicVolumeCropper.",
            "The TOADS equality check confirms saved crops copy the already-transferred full-grid maps; it does not validate the earlier oblique-grid transfer against MIPAV ReorientVolume.",
            "This comparison does not validate cortical separation or whole-workflow equivalence.",
        ],
    }
    report_path = out / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    summary = {
        "report": str(report_path),
        "bounds_inclusive": report["bounds_inclusive"],
        "crop_shape_xyz": report["crop_shape_xyz"],
        "same_grid_toads_identity_checks": same_grid_rows,
        "crop_outputs": rows,
        "all_crop_outputs_exact": report["all_crop_outputs_exact"],
    }
    print(json.dumps(summary, indent=2))
    return 0 if report["all_crop_outputs_exact"] and all(
        bool(row["exact_voxel_copy_from_same_grid_input"])
        for row in same_grid_rows.values()
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
