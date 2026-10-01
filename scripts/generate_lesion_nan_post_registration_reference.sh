#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture="$repo_root/tests/data/synthetic/lesion_post_registration_nan"
reference="$repo_root/tests/java_reference/LesionNaNPostRegistrationBranchReferenceHarness.java"
classes="$repo_root/build/java-reference/lesion-post-registration-nan"
mkdir -p "$fixture" "$classes"
javac -d "$classes" "$reference"
java -cp "$classes" LesionNaNPostRegistrationBranchReferenceHarness > "$fixture/java_reference.json"
"$repo_root/.venv/bin/python" - "$fixture/java_reference.json" "$reference" <<'PY'
import hashlib
import json
import pathlib
import sys

out, harness = map(pathlib.Path, sys.argv[1:])
reference = json.loads(out.read_text())
reference["harness_sha256"] = hashlib.sha256(harness.read_bytes()).hexdigest()
out.write_text(json.dumps(reference, separators=(",", ":")) + "\n")
PY
