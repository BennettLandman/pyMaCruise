#!/usr/bin/env python3
"""Regenerate selected-Java bits for the exact P39 initial-thinning state."""

from __future__ import annotations

import hashlib
import json
import struct
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/lesion_initial_thinning_composition"
OUTPUT = ROOT / "tests/data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json"
HARNESS = ROOT / "tests/java_reference/LesionPostThinningStatisticsReferenceHarness.java"
POWER = ROOT / "tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"


def _unsigned_json(value: object) -> object:
    if isinstance(value, list):
        return [_unsigned_json(item) for item in value]
    if isinstance(value, int):
        return value & 0xFFFFFFFF
    return value


def _input_bytes() -> bytes:
    spec = json.loads((FIXTURE / "fixture.json").read_text(encoding="utf-8"))
    reference = json.loads((FIXTURE / "java_reference.json").read_text(encoding="utf-8"))
    frames = np.asarray(nib.load(str(FIXTURE / "initial_thinning_inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    if frames.shape != (17, 17, 17, 64):
        raise ValueError(f"unexpected source phantom shape {frames.shape}")
    memberships = np.asarray(reference["class_membership_bits"], dtype=np.uint32).view(np.float32).reshape(17, 17, 17, 13)
    segmentation = np.asarray(reference["final_segmentation"], dtype=np.int32).reshape(17, 17, 17)
    chunks = [np.asarray(spec["class_centroids"], dtype=">f4").tobytes()]
    chunks.append(np.asarray(spec["lesion_centroid"], dtype=">f4").tobytes())
    # This synthetic profile has no BLACKHOLE class; zero is the explicit
    # fixture initializer and the selected method sets an empty estimate to 0.
    chunks.append(np.zeros(2, dtype=">f4").tobytes())
    for x in range(17):
        for y in range(17):
            for z in range(17):
                chunks.append(np.asarray(frames[x, y, z, 0:2], dtype=">f4").tobytes())
                chunks.append(np.asarray(memberships[x, y, z, :], dtype=">f4").tobytes())
                chunks.append(struct.pack(">i", int(segmentation[x, y, z])))
                chunks.append(struct.pack(">f", float(frames[x, y, z, 60])))
                chunks.append(struct.pack(">f", float(frames[x, y, z, 61])))
    return b"".join(chunks)


def main() -> None:
    input_bytes = _input_bytes()
    with tempfile.TemporaryDirectory(prefix="pymacruise-p39-stats-") as temp_name:
        temp = Path(temp_name)
        classes = temp / "classes"
        classes.mkdir()
        input_path = temp / "p39-post-thinning.bin"
        input_path.write_bytes(input_bytes)
        subprocess.run(
            ["javac", "-d", str(classes), str(POWER), str(HARNESS)],
            check=True,
        )
        completed = subprocess.run(
            ["java", "-cp", str(classes), "LesionPostThinningStatisticsReferenceHarness", str(input_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    result = json.loads(completed.stdout)
    for key in ("class", "aux_centroids", "class_std", "lesion_std", "blackhole_std", "weights"):
        result[key] = _unsigned_json(result[key])
    result["harness_sha256"] = hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    result["power_table_file_sha256"] = hashlib.sha256(POWER.read_bytes()).hexdigest()
    result["input_sha256"] = hashlib.sha256(input_bytes).hexdigest()
    result["upstream_membership_reference_sha256"] = hashlib.sha256(
        (FIXTURE / "java_reference.json").read_bytes()
    ).hexdigest()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
