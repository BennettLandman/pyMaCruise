#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT
javac -d "$build_dir" "$repo_root/tests/java_reference/LesionToadsStandard13VoxelReferenceHarness.java"
java -cp "$build_dir" LesionToadsStandard13VoxelReferenceHarness \
  > "$repo_root/tests/data/synthetic/lesion_standard13_voxel/java_reference.json"
python3 - "$build_dir" "$repo_root" <<'PY'
import hashlib
import json
import pathlib
import subprocess
import sys

build_dir = pathlib.Path(sys.argv[1])
repo_root = pathlib.Path(sys.argv[2])
harness = repo_root / "tests/java_reference/LesionToadsStandard13VoxelReferenceHarness.java"
source_sha = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
cases = {
    "valid_gm": (-1, 2, 10, True),
    "invalid_gm": (-1, -1, 10, True),
    "zero_gm": (-1, 0, 10, True),
    "zero_gm_without_brainstem": (-1, 0, -1, True),
    "ventricle_precedence": (1, 2, 10, True),
    "ventricle_only": (1, -1, 10, True),
    "gm_without_brainstem": (-1, 2, -1, True),
}
references = {}
for name, args in cases.items():
    output = subprocess.check_output(
        ["java", "-cp", str(build_dir), "LesionToadsStandard13VoxelReferenceHarness", *map(str, args)],
        text=True,
    )
    references[name] = json.loads(output)
payload = {
    "harness_sha256": hashlib.sha256(harness.read_bytes()).hexdigest(),
    "source_sha256": source_sha,
    "cases": references,
}
path = repo_root / "tests/data/synthetic/lesion_standard13_voxel/gm_flag_java_reference.json"
path.write_text(json.dumps(payload, indent=2) + "\n")
PY
