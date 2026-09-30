#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
source_file="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
wrapper_file="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmLesionToads.java"
source_sha="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
wrapper_sha="556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
test_source_sha=$(shasum -a 256 "$source_file" | awk '{print $1}')
test_wrapper_sha=$(shasum -a 256 "$wrapper_file" | awk '{print $1}')
if [ "$test_source_sha" != "$source_sha" ] || [ "$test_wrapper_sha" != "$wrapper_sha" ]; then
    echo "Selected LesionTOADS source or wrapper hash changed" >&2
    exit 1
fi
fixture="$repo_root/tests/data/synthetic/lesion_statistics_membership_bridge"
"$repo_root/.venv/bin/python" "$fixture/generate_fixture.py"
classes="$repo_root/build/java-reference/lesion-statistics-membership-bridge"
mkdir -p "$classes"
javac -d "$classes" "$repo_root/tests/java_reference/LesionStatisticsMembershipBridgeReferenceHarness.java"
java -cp "$classes" LesionStatisticsMembershipBridgeReferenceHarness > "$fixture/java_reference.json"
python - "$fixture/java_reference.json" "$repo_root/tests/java_reference/LesionStatisticsMembershipBridgeReferenceHarness.java" <<'PY'
import hashlib
import json
import pathlib
import sys

destination = pathlib.Path(sys.argv[1])
reference = json.loads(destination.read_text())
reference["source_sha256"] = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
reference["wrapper_sha256"] = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
reference["harness_sha256"] = hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest()
destination.write_text(json.dumps(reference, separators=(",", ":")) + "\n")
PY
