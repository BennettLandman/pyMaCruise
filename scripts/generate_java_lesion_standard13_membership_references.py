#!/usr/bin/env python3
"""Regenerate the fixed-input selected-Java standard-13 membership oracles."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
JAVA = ROOT / "tests" / "java_reference"
POWER_TABLE = JAVA / "upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
HARNESS = JAVA / "LesionToadsStandard13VoxelReferenceHarness.java"
VOLUME_HARNESS = JAVA / "LesionToadsStandard13VolumeReferenceHarness.java"
CLASSES = ROOT / "build/java-reference/lesion-standard13-classes"
VOXEL_REF = ROOT / "tests/data/synthetic/lesion_standard13_voxel/java_reference.json"
GM_REF = ROOT / "tests/data/synthetic/lesion_standard13_voxel/gm_flag_java_reference.json"
VOLUME_REF = ROOT / "tests/data/synthetic/lesion_standard13_volume/java_reference.json"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def main() -> None:
    CLASSES.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["javac", "-d", str(CLASSES), str(POWER_TABLE), str(HARNESS), str(VOLUME_HARNESS)],
        check=True,
    )
    harness_sha = hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    cases = {
        "valid_gm": ("-1", "2", "10", "true"),
        "invalid_gm": ("-1", "-1", "10", "true"),
        "zero_gm": ("-1", "0", "10", "true"),
        "zero_gm_without_brainstem": ("-1", "0", "-1", "true"),
        "ventricle_precedence": ("1", "2", "10", "true"),
        "ventricle_only": ("1", "-1", "10", "true"),
        "gm_without_brainstem": ("-1", "2", "-1", "true"),
    }
    outputs = {
        name: json.loads(
            subprocess.check_output(
                ["java", "-cp", str(CLASSES), "LesionToadsStandard13VoxelReferenceHarness", *args],
                text=True,
            )
        )
        for name, args in cases.items()
    }
    VOXEL_REF.write_text(
        subprocess.check_output(
            ["java", "-cp", str(CLASSES), "LesionToadsStandard13VoxelReferenceHarness"],
            text=True,
        )
    )
    VOLUME_REF.write_text(
        subprocess.check_output(
            ["java", "-cp", str(CLASSES), "LesionToadsStandard13VolumeReferenceHarness"],
            text=True,
        )
    )
    GM_REF.write_text(
        json.dumps(
            {"harness_sha256": harness_sha, "source_sha256": SOURCE_SHA, "cases": outputs},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
