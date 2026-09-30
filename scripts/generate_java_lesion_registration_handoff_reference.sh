#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture="$repo_root/tests/data/synthetic/lesion_registration_handoff"
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
classes="$repo_root/build/java-reference/lesion-registration-handoff"

"$repo_root/.venv/bin/python" "$fixture/generate_fixture.py"
mkdir -p "$classes"
javac -d "$classes" \
  "$reference_root/stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$upstream/utilities/Numerics.java" \
  "$upstream/structures/RotationMatrix.java" \
  "$upstream/methods/ParametricTransform.java" \
  "$upstream/libraries/ImageFunctions.java" \
  "$reference_root/LesionRegistrationHandoffReferenceHarness.java"
java -cp "$classes" LesionRegistrationHandoffReferenceHarness \
  "$fixture/no_field_memberships.txt" "$fixture/with_field_memberships.txt" \
  > "$fixture/java_reference.json"
python3 - "$fixture/java_reference.json" "$reference_root/LesionRegistrationHandoffReferenceHarness.java" <<'PY'
import hashlib
import json
import pathlib
import sys

output = pathlib.Path(sys.argv[1])
reference = json.loads(output.read_text())
reference["atlas_source_sha256"] = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
reference["wrapper_source_sha256"] = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
reference["harness_sha256"] = hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest()
output.write_text(json.dumps(reference, separators=(",", ":")) + "\n")
PY
