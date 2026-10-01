#!/usr/bin/env python3
"""Generate the selected-Java reference for the post-thinning GM update."""

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
OUTPUT = ROOT / "tests/data/synthetic/lesion_post_thinning_gm_membership_update/java_reference.json"
HARNESS = ROOT / "tests/java_reference/LesionPostThinningStandard13MembershipReferenceHarness.java"
POWER = ROOT / "tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"


def _f32(bits: object) -> np.ndarray:
    return np.asarray(bits, dtype=np.uint32).view(np.float32)


def _input_bytes() -> bytes:
    spec = json.loads((FIXTURE / "fixture.json").read_text(encoding="utf-8"))
    first = json.loads((FIXTURE / "java_reference.json").read_text(encoding="utf-8"))
    statistics = json.loads(
        (ROOT / "tests/data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_text(encoding="utf-8")
    )
    if first["source_sha256"] != statistics["source_sha256"]:
        raise ValueError("first-pass and statistics Java references use different selected sources")
    if statistics["upstream_membership_reference_sha256"] != hashlib.sha256((FIXTURE / "java_reference.json").read_bytes()).hexdigest():
        raise ValueError("statistics reference is not pinned to the current first-pass membership reference")
    frames = np.asarray(nib.load(str(FIXTURE / "initial_thinning_inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    if frames.shape != (17, 17, 17, 64) or nib.aff2axcodes(nib.load(str(FIXTURE / "initial_thinning_inputs_ras.nii.gz")).affine) != ("R", "A", "S"):
        raise ValueError("the selected P39 RAS fixture changed shape or orientation")
    membership = _f32(first["class_membership_bits"]).reshape(17, 17, 17, 13)
    lesion = _f32(first["lesion_membership_bits"]).reshape(17, 17, 17)
    black = _f32(first["black_hole_membership_bits"]).reshape(17, 17, 17)
    outlier = _f32(first["outlier_membership_bits"]).reshape(17, 17, 17)
    relation = _f32(first["lesion_relation_bits"]).reshape(17, 17, 17)
    chunks: list[bytes] = []
    for key in ("class",):
        chunks.append(np.asarray(_f32(statistics[key]), dtype=">f4").reshape(2, 13).tobytes())
    chunks.append(np.asarray(_f32(statistics["class_std"]), dtype=">f4").reshape(2, 13).tobytes())
    # The standard-13 voxel adapter consumes the lesion centroid only; the
    # fixture has no BLACKHOLE class and keeps that initialized map at zero.
    chunks.append(np.asarray(_f32(statistics["aux_centroids"])[0], dtype=">f4").tobytes())
    chunks.append(np.asarray(_f32(statistics["lesion_std"]), dtype=">f4").tobytes())
    chunks.append(np.asarray(_f32(statistics["weights"]), dtype=">f4").tobytes())
    chunks.append(np.asarray(spec["shape_priors"], dtype=">f4").tobytes())
    chunks.append(np.asarray(spec["smoothing_priors"], dtype=">f4").tobytes())
    for x in range(17):
        for y in range(17):
            for z in range(17):
                chunks.append(np.asarray(frames[x, y, z, 0:2], dtype=">f4").tobytes())
                chunks.append(np.asarray(frames[x, y, z, 2:15], dtype=">f4").tobytes())
                chunks.append(np.asarray(frames[x, y, z, 15:28] > 0.5, dtype=">i4").tobytes())
                chunks.append(np.asarray(membership[x, y, z], dtype=">f4").tobytes())
                chunks.append(struct.pack(">4f", float(lesion[x, y, z]), float(black[x, y, z]), float(outlier[x, y, z]), float(relation[x, y, z])))
                chunks.append(struct.pack(">3i", int(frames[x, y, z, 44]), int(frames[x, y, z, 45]), int(frames[x, y, z, 46])))
    return b"".join(chunks)


def main() -> None:
    FIXTURE_OUT = OUTPUT.parent
    FIXTURE_OUT.mkdir(parents=True, exist_ok=True)
    statistics = json.loads(
        (ROOT / "tests/data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_text(encoding="utf-8")
    )
    payload = _input_bytes()
    with tempfile.TemporaryDirectory(prefix="pymacruise-p39-gm-membership-") as temp_name:
        temp = Path(temp_name)
        classes = temp / "classes"
        classes.mkdir()
        input_path = temp / "state.bin"
        input_path.write_bytes(payload)
        subprocess.run(["javac", "-d", str(classes), str(POWER), str(HARNESS)], check=True)
        result = subprocess.run(
            ["java", "-cp", str(classes), "LesionPostThinningStandard13MembershipReferenceHarness", str(input_path)],
            check=True, capture_output=True, text=True,
        )
    reference = json.loads(result.stdout)
    reference["wrapper_sha256"] = statistics["wrapper_sha256"]
    reference["power_table_file_sha256"] = statistics["power_table_file_sha256"]
    reference["source_call_order"] = [
        "propagateThinning:heapLoop",
        "computeCentroids",
        "computeVariances",
        "computeChannelWeights",
        "computeMemberships(GMflag=true)",
    ]
    reference["harness_sha256"] = hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    reference["input_sha256"] = hashlib.sha256(payload).hexdigest()
    reference["upstream_initial_membership_reference_sha256"] = hashlib.sha256((FIXTURE / "java_reference.json").read_bytes()).hexdigest()
    reference["upstream_statistics_reference_sha256"] = hashlib.sha256((ROOT / "tests/data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_bytes()).hexdigest()
    reference["fixture_sha256"] = hashlib.sha256((FIXTURE / "initial_thinning_inputs_ras.nii.gz").read_bytes()).hexdigest()
    OUTPUT.write_text(json.dumps(reference, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
