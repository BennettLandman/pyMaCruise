#!/usr/bin/env python3
"""Compare saved PLUS Dura transfer with selected JIST UNCHANGED behavior.

This replays ``ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR)`` on the
saved Dura channels and T1 reference, then compares against Python's explicit
``jist-unchanged`` transfer mode. It is a source-compatibility check; the
recommended physical-affine transfer is a separate behavior.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from compare_human_jist_toads_transfer import ensure_compiled  # noqa: E402
from pymacruise.spatial_transfer import assign_reference_geometry_without_resampling  # noqa: E402


HARNESS = REPO / "tests/java_reference/ReorientVolumeHumanTransferReferenceHarness.java"
CHANNELS = ("csf", "gm", "wm_fill", "wm_mask")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare(workflow_root: Path, output_dir: Path) -> dict:
    workflow_root = workflow_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    transfer = json.loads((workflow_root / "native-dura/native-dura-transfer.json").read_text())
    if transfer.get("transfer_mode") != "physical-affine":
        raise ValueError("expected a workflow whose recommended physical-affine outputs are saved")
    reference_path = Path(transfer["inputs"]["t1_reference"]["path"])
    reference = nib.load(str(reference_path))
    if nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("T1 reference must be RAS-oriented")
    classpath = ensure_compiled()
    java_home = output_dir / "java-home"
    java_home.mkdir(exist_ok=True)
    results = {}
    for channel in CHANNELS:
        source_path = Path(transfer["inputs"][f"dura_{channel}"]["path"])
        if sha256(source_path) != transfer["inputs"][f"dura_{channel}"]["sha256"]:
            raise ValueError(f"saved Dura {channel} input hash mismatch")
        source = nib.load(str(source_path))
        if nib.aff2axcodes(source.affine) != ("R", "A", "S"):
            raise ValueError(f"Dura {channel} input must be RAS-oriented")
        base = output_dir / channel
        command = [
            "java", "-Djava.awt.headless=true", f"-Duser.home={java_home}",
            "-Xms128m", "-Xmx2g", "-cp", classpath,
            "ReorientVolumeHumanTransferReferenceHarness", str(reference_path),
            channel, str(source_path), f"{base}.java.f32be,{base}.java.json",
        ]
        with (output_dir / "java-harness.log").open("a") as log:
            subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
        geometry = json.loads(Path(f"{base}.java.json").read_text())
        expected_image = assign_reference_geometry_without_resampling(source, reference)
        expected = np.asarray(expected_image.dataobj, dtype=np.float32)
        if tuple(geometry["output_model_geometry_after_reference_header"]["axis_orientation"]) != (2, 4, 5):
            raise ValueError(f"unexpected MIPAV orientation on {channel}: {geometry}")
        raw = np.fromfile(f"{base}.java.f32be", dtype=">f4")
        if raw.size != expected.size:
            raise ValueError(f"Java returned wrong voxel count for {channel}: {raw.size} vs {expected.size}")
        actual = np.flip(raw.reshape(expected.shape, order="F"), axis=1)
        delta = np.abs(actual.astype(np.float64) - expected.astype(np.float64))
        results[channel] = {
            "source_shape": list(source.shape),
            "reference_shape": list(reference.shape),
            "java_extents": geometry["output_model_geometry_after_reference_header"]["extents"],
            "python_shape": list(expected.shape),
            "bitwise_equal_float32": bool(np.array_equal(actual, expected)),
            "changed_voxels": int(np.count_nonzero(actual != expected)),
            "max_absolute_error": float(delta.max(initial=0.0)),
            "mean_absolute_error": float(delta.mean()),
        }
    return {
        "case": "corrected atlas-aligned synthetic PLUS Dura transfer",
        "synthetic_only": True,
        "source_method": "JIST ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR) then reference-header assignment",
        "python_method": "assign_reference_geometry_without_resampling (jist-unchanged)",
        "reference": {"path": str(reference_path), "sha256": sha256(reference_path), "shape": list(reference.shape)},
        "coordinate_note": "The source method preserves the Dura sample extents and values, then assigns reference geometry; it does not physically resample into reference dimensions.",
        "results": results,
        "all_channels_exact": all(item["bitwise_equal_float32"] for item in results.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.workflow_root, args.output_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_channels_exact"]:
        raise SystemExit("Python legacy transfer differs from selected Java JIST")


if __name__ == "__main__":
    main()
