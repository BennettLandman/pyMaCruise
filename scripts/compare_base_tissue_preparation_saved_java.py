#!/usr/bin/env python3
"""Compare saved base MaSeg/TOADS fusion outputs with selected Java source.

The preparation directory must contain ``provenance.json`` written by the
``base-tissue-preparation`` command. Six cropped images are held fixed and
sent to the existing headless ``CombineMAandTOADS`` source harness. This
compares fusion only; it does not execute or validate transfer/cropping.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np


INPUT_KEYS = (
    "atlas_gm_cropped",
    "atlas_wm_cropped",
    "atlas_brain_tissue_mask_cropped",
    "toads_gm_cropped",
    "toads_wmfill_cropped",
    "toads_mask_cropped",
)
OUTPUT_KEYS = {"gm": "fused_gm", "wm": "fused_wmfill", "mask": "fused_mask"}
EXPECTED_ORDER = [
    "corticalSeparationFromMaSeg",
    "nearestNeighborTransferToMaSegReference",
    "sharedCropFromPositiveMaSegGM",
    "combineMaSegAtlasAndTOADS",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
    if parameters.get("stage_order") != EXPECTED_ORDER:
        raise SystemExit("Preparation stage order does not match the selected V3 route")

    images = [nib.load(provenance["outputs"][key]) for key in INPUT_KEYS]
    shape = tuple(int(size) for size in images[0].shape)
    if len(shape) != 3 or any(tuple(image.shape) != shape for image in images):
        raise SystemExit("Held-fixed fusion inputs do not share one 3D shape")
    reference_affine = images[0].affine
    if any(
        not np.allclose(image.affine, reference_affine, rtol=0, atol=1e-5)
        for image in images[1:]
    ):
        raise SystemExit("Held-fixed fusion inputs do not share one physical grid")
    arrays = [np.asarray(image.dataobj, dtype=np.float32) for image in images]
    binary_path = out / "tissue-arrays.bin"
    with binary_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for array in arrays:
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="F"))

    for name, key in OUTPUT_KEYS.items():
        saved = nib.load(provenance["outputs"][key])
        if tuple(saved.shape) != shape or not np.allclose(
            saved.affine, reference_affine, rtol=0, atol=1e-5
        ):
            raise SystemExit(f"Saved Python fusion output {key} is off the input grid")
        np.save(out / f"python-{name}.npy", np.asarray(saved.dataobj))

    source_root = repo / "tests/java_reference/upstream"
    stubs_root = repo / "tests/java_reference/stubs"
    harness = repo / "tests/java_reference/plus_tissue_update/PlusTissueUpdateSourceHarness.java"
    classes = out / "classes"
    classes.mkdir(exist_ok=True)
    javac = shutil.which("javac")
    java = shutil.which("java")
    if not javac or not java:
        raise SystemExit("javac and java are required for the preserved Java comparison")
    subprocess.run(
        [
            javac,
            "-sourcepath",
            f"{stubs_root}:{source_root}",
            "-d",
            str(classes),
            str(source_root / "edu/vanderbilt/masi/plugins/CRUISE/utilities/MASImorph.java"),
            str(source_root / "edu/vanderbilt/masi/plugins/CRUISE/MaCRUISE/CombineMAandTOADS.java"),
            str(harness),
        ],
        check=True,
    )
    watchdog = repo / "scripts/run_with_rss_watchdog.py"
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
            "edu.vanderbilt.masi.plugins.CRUISE.MaCRUISE.PlusTissueUpdateSourceHarness",
            str(binary_path),
            str(out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    rows: dict[str, dict[str, object]] = {}
    for name in OUTPUT_KEYS:
        with gzip.open(out / f"java-{name}.raw.gz", "rb") as stream:
            java_shape = struct.unpack(">III", stream.read(12))
            java_values = np.frombuffer(stream.read(), dtype=">f4").astype(np.float32)
        if java_shape != shape or java_values.size != int(np.prod(shape)):
            raise SystemExit(f"Java {name} output has an unexpected shape")
        java_array = java_values.reshape(shape, order="C")
        python_array = np.load(out / f"python-{name}.npy").astype(np.float32, copy=False)
        delta = np.abs(java_array.astype(np.float64) - python_array.astype(np.float64))
        rows[name] = {
            "exact": bool(np.array_equal(java_array, python_array)),
            "changed_voxels": int(np.count_nonzero(java_array != python_array)),
            "voxel_count": int(java_array.size),
            "maximum_absolute_error": float(delta.max(initial=0)),
            "mean_absolute_error": float(delta.mean()),
        }

    report = {
        "comparison": "selected source Java CombineMAandTOADS.execute() versus saved Python V3 base-preparation fusion outputs, with six cropped inputs held fixed",
        "source_files": [
            str(source_root / "edu/vanderbilt/masi/plugins/CRUISE/MaCRUISE/CombineMAandTOADS.java"),
            str(source_root / "edu/vanderbilt/masi/plugins/CRUISE/utilities/MASImorph.java"),
        ],
        "source_settings": {"erosion": "3D center plus six neighbors", "gm_confident_level": 0.001},
        "java_harness": str(harness),
        "java_resource_limit": "-Xmx768m with 2 GiB sampled-RSS watchdog",
        "watchdog_output": java_result.stdout.strip(),
        "preparation_provenance": str(provenance_path),
        "preparation_provenance_sha256": _sha256(provenance_path),
        "remove_cerebellum": False,
        "preparation_stage_order": EXPECTED_ORDER,
        "input_files": {
            key: {
                "path": provenance["outputs"][key],
                "sha256": _sha256(Path(provenance["outputs"][key])),
            }
            for key in INPUT_KEYS
        },
        "input_binary_sha256": _sha256(binary_path),
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "outputs": rows,
        "all_outputs_exact": all(bool(row["exact"]) for row in rows.values()),
        "scope_limitations": [
            "This comparison validates fusion only on held-fixed Python-prepared cropped inputs.",
            "Cortical separation, Java ReorientVolume transfer, CropBasedOnGM/CubicVolumeCropper bounds and voxel copying, and full-workflow equivalence remain unvalidated by this comparison.",
        ],
    }
    report_path = out / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"report": str(report_path), **{key: report[key] for key in ("shape_xyz", "outputs", "all_outputs_exact")}}, indent=2))
    return 0 if report["all_outputs_exact"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
