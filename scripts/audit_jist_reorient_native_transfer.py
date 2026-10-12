#!/usr/bin/env python3
"""Audit JIST ReorientVolume versus physical-RAS native-grid transfer.

This source audit reports output extents and whether JIST preserves the input
voxel array when source and reference are already RAS with matching axis
orientation. It also records the Python physical-affine adapter's full
reference-grid geometry. It does not claim the two operators are equivalent
for differing grids.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from compare_human_jist_toads_transfer import (
    BUILD,
    HARNESS,
    REORIENT,
    classpath,
    ensure_compiled,
    image_geometry,
)
from pymacruise.spatial_transfer import transfer_to_reference_space


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_map(value: str) -> tuple[str, Path]:
    label, separator, raw_path = value.partition("=")
    if not separator or not label or not raw_path:
        raise argparse.ArgumentTypeError("maps must use LABEL=IMAGE.nii.gz")
    path = Path(raw_path).expanduser().resolve()
    return label, path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="native reference RAS .nii.gz")
    parser.add_argument("--map", action="append", required=True, type=parse_map, dest="maps", help="source RAS map, LABEL=IMAGE.nii.gz; repeatable")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    reference_path = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not reference_path.is_file() or reference_path.is_symlink():
        raise FileNotFoundError(f"Reference is not a regular file: {reference_path}")
    if len({label for label, _ in args.maps}) != len(args.maps):
        raise ValueError("map labels must be unique")
    reference = nib.load(str(reference_path))
    if len(reference.shape) != 3 or nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("reference must be a 3D RAS NIfTI")
    for _, source_path in args.maps:
        if not source_path.is_file() or source_path.is_symlink():
            raise FileNotFoundError(f"Source map is not a regular file: {source_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    java_home = output_dir / "java-home"
    java_home.mkdir(parents=True, exist_ok=True)
    cp = ensure_compiled()
    java_args = [
        "java", "-Djava.awt.headless=true", f"-Duser.home={java_home}",
        "-Xms128m", "-Xmx2g", "-cp", cp,
        "ReorientVolumeHumanTransferReferenceHarness", str(reference_path),
    ]
    for label, source_path in args.maps:
        java_args.extend((label, str(source_path), f"{output_dir / label}.java.f32be,{output_dir / label}.java.json"))
    with (output_dir / "java-harness.log").open("w", encoding="utf-8") as log:
        subprocess.run(java_args, cwd=Path(__file__).resolve().parents[1], stdout=log, stderr=subprocess.STDOUT, check=True)

    records = {}
    for label, source_path in args.maps:
        source = nib.load(str(source_path))
        if len(source.shape) != 3 or nib.aff2axcodes(source.affine) != ("R", "A", "S"):
            raise ValueError(f"{label} must be a 3D RAS NIfTI")
        java_metadata_path = output_dir / f"{label}.java.json"
        java_metadata = json.loads(java_metadata_path.read_text(encoding="utf-8"))
        output_geometry = java_metadata["output_model_geometry_after_reference_header"]
        java_extents = tuple(int(value) for value in output_geometry["extents"])
        java_linear = np.fromfile(output_dir / f"{label}.java.f32be", dtype=">f4")
        if java_linear.size != int(np.prod(java_extents)):
            raise ValueError(f"{label}: Java raw output does not match its recorded extents")
        # MIPAV's A2P storage axis runs opposite the NIfTI +A array order.
        java_ras = np.flip(java_linear.reshape(java_extents, order="F"), axis=1)
        source_values = np.asarray(source.dataobj, dtype=np.float32)
        same_array = source.shape == java_extents and np.array_equal(java_ras, source_values)
        python_physical = transfer_to_reference_space(source, reference, cval=0.0)
        records[label] = {
            "source": {"path": str(source_path), "sha256": sha256(source_path), "geometry": image_geometry(source)},
            "java_output_geometry_after_reference_header": output_geometry,
            "java_output_extents_xyz": list(java_extents),
            "java_voxel_count": int(java_linear.size),
            "java_array_matches_source_voxels_exactly": bool(same_array),
            "java_data_sha256": sha256(output_dir / f"{label}.java.f32be"),
            "python_physical_transfer_geometry": image_geometry(python_physical),
            "output_shape_matches": list(java_extents) == list(python_physical.shape),
        }

    report = {
        "case": "selected JIST ReorientVolume versus PyMaCRUISE physical-affine transfer",
        "evidence_category": "held-fixed source geometry audit",
        "java_source": str(REORIENT),
        "java_source_sha256": sha256(REORIENT),
        "java_harness_sha256": sha256(HARNESS),
        "reference": {"path": str(reference_path), "sha256": sha256(reference_path), "geometry": image_geometry(reference)},
        "maps": records,
        "interpretation": "JIST uses Resolution.UNCHANGED and may preserve source extents; Python samples on the full reference grid. Different output extents are reported, not silently compared voxelwise.",
    }
    report_path = output_dir / "geometry-audit.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report={report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
