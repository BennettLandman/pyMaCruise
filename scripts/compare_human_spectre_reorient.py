#!/usr/bin/env python3
"""Replay source JIST axial reorientation on a human N4 image.

This is a diagnostic geometry check for historical SPECTRE preprocessing. It
does not run skull stripping, atlas registration, or the selected native-RAS
PyMaCRUISE workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path("/Users/landmaba/dev/mipav/sample-macruise/output")
DEFAULT_N4 = DEFAULT_SOURCE / "ToadsSeperate/N4Correct/T1_N4Corrected.nii.gz"
DEFAULT_REFERENCE = DEFAULT_SOURCE / "Raw/T1.nii.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-n4", type=Path, default=DEFAULT_N4)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "build/human-spectre-reorient-2026-10-12",
    )
    args = parser.parse_args()
    source_n4 = args.source_n4.resolve()
    reference = args.reference.resolve()
    output_dir = args.output_dir.resolve()
    for path in (source_n4, reference):
        if not path.is_file():
            parser.error(f"input does not exist: {path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Reuse the already-audited headless Java/MIPAV classpath setup. Build a
    # local copy of the existing harness with only the selected SPECTRE call
    # changed to AXIAL/UNCHANGED/WINDOWED_SINC.
    sys.path.insert(0, str(ROOT / "scripts"))
    import compare_human_jist_toads_transfer as transfer_audit  # noqa: PLC0415

    transfer_audit.ensure_compiled()
    harness_text = transfer_audit.HARNESS.read_text(encoding="utf-8")
    old_call = (
        "ReorientVolume.solve(\n"
        "                    source, reference, Resolution.UNCHANGED, Interpolation.NEAREST_NEIGHBOR)"
    )
    new_call = (
        "ReorientVolume.solve(\n"
        "                    source, reference, Orientation.AXIAL, Resolution.UNCHANGED, "
        "Interpolation.WINDOWED_SINC)"
    )
    if harness_text.count(old_call) != 1:
        raise RuntimeError("Reference harness call changed; refusing an implicit rewrite")
    harness_text = harness_text.replace(
        "import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume.Resolution;",
        "import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume.Resolution;\n"
        "import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume.Orientation;",
    )
    harness_text = harness_text.replace(
        "ReorientVolumeHumanTransferReferenceHarness",
        "ReorientVolumeAxialReferenceHarness",
    ).replace(old_call, new_call)
    harness_path = output_dir / "ReorientVolumeAxialReferenceHarness.java"
    harness_path.write_text(harness_text, encoding="utf-8")
    classpath = transfer_audit.classpath()
    subprocess.run(
        ["javac", "-classpath", classpath, "-d", str(transfer_audit.BUILD / "classes"), str(harness_path)],
        cwd=ROOT,
        check=True,
    )

    raw_path = output_dir / "java-reoriented-n4.f32be"
    geometry_path = output_dir / "java-reoriented-n4.json"
    java_home = output_dir / "java-home"
    java_home.mkdir(exist_ok=True)
    command = [
        "java", "-Djava.awt.headless=true", f"-Duser.home={java_home}",
        "-Xms128m", "-Xmx2g", "-cp", classpath,
        "ReorientVolumeAxialReferenceHarness", str(reference), "n4",
        str(source_n4), f"{raw_path},{geometry_path}",
    ]
    java_log = output_dir / "java-harness.log"
    with java_log.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        tail = "\n".join(java_log.read_text(errors="replace").splitlines()[-35:])
        raise RuntimeError(f"JIST axial reorientation harness failed:\n{tail}")

    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    output_geometry = geometry["output_model_geometry_before_reference_header"]
    extents = tuple(int(value) for value in output_geometry["extents"])
    java_values = np.fromfile(raw_path, dtype=">f4")
    if java_values.size != int(np.prod(extents)):
        raise ValueError("Java output voxel count does not match its geometry")
    java_values = java_values.reshape(extents, order="F")
    source_image = nib.load(str(source_n4))
    if tuple(source_image.shape) != extents:
        raise ValueError("Java axial reorientation changed the source N4 shape")
    source_values = np.asarray(source_image.dataobj, dtype=np.float32)
    # Source N4 is RAS NIfTI; MIPAV's R2L/A2P/I2S axial buffer reverses x/y.
    expected_values = np.flip(source_values, axis=(0, 1))
    delta = np.abs(java_values.astype(np.float64) - expected_values.astype(np.float64))
    correlation = float(np.corrcoef(java_values.ravel(), expected_values.ravel())[0, 1])
    report = {
        "scope": "selected JIST ReorientVolume axial orientation replay; not SPECTRE mask or atlas-registration validation",
        "source_n4": {"path": str(source_n4), "sha256": sha256(source_n4)},
        "reference": {"path": str(reference), "sha256": sha256(reference)},
        "java_source": str(transfer_audit.REORIENT),
        "java_source_sha256": sha256(transfer_audit.REORIENT),
        "operation": {"orientation": "AXIAL", "resolution": "UNCHANGED", "interpolation": "WINDOWED_SINC"},
        "java_output_geometry": output_geometry,
        "expected_ras_to_lps_array_conversion": "flip axes 0 and 1",
        "comparison": {
            "shape": list(extents),
            "changed_voxels": int(np.count_nonzero(java_values != expected_values)),
            "mean_absolute_error": float(delta.mean()),
            "maximum_absolute_error": float(delta.max(initial=0.0)),
            "correlation": correlation,
        },
        "java_log": str(java_log),
    }
    report_path = output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), **report["comparison"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
