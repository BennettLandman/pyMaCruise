#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture="$repo_root/tests/data/synthetic/lesion_post_thinning_registration"
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
classes="$repo_root/build/java-reference/lesion-post-thinning-registration"
mkdir -p "$classes"

"$repo_root/.venv/bin/python" "$fixture/generate_fixture.py"
javac -d "$classes" \
  "$reference_root/stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$upstream/utilities/Numerics.java" \
  "$upstream/structures/RotationMatrix.java" \
  "$upstream/methods/ParametricTransform.java" \
  "$upstream/libraries/ImageFunctions.java" \
  "$reference_root/LesionPostThinningRegistrationReferenceHarness.java"
java -cp "$classes" LesionPostThinningRegistrationReferenceHarness \
  "$fixture/registration_inputs_ras.nii.gz" > "$fixture/java_reference.json"
"$repo_root/.venv/bin/python" - "$fixture" "$reference_root/LesionPostThinningRegistrationReferenceHarness.java" \
  "$upstream/methods/ParametricTransform.java" "$upstream/libraries/ImageFunctions.java" \
  "$upstream/structures/RotationMatrix.java" <<'PY'
import hashlib
import json
import pathlib
import sys

fixture = pathlib.Path(sys.argv[1])
reference = json.loads((fixture / "java_reference.json").read_text())
reference["algorithm_source_sha256"] = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
reference["atlas_source_sha256"] = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
for name, source in zip(("harness", "transform", "interpolation", "rotation"), sys.argv[2:]):
    reference[f"{name}_sha256"] = hashlib.sha256(pathlib.Path(source).read_bytes()).hexdigest()
(fixture / "java_reference.json").write_text(json.dumps(reference, separators=(",", ":")) + "\n")
PY
