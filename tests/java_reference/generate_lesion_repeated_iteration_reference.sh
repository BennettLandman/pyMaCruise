#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
HARNESS="$ROOT/tests/java_reference/LesionToadsRepeatedIterationReferenceHarness.java"
OUT="$ROOT/tests/data/synthetic/lesion_repeated_iteration/java_reference.json"
LOOP_STATES="$ROOT/tests/data/synthetic/lesion_repeated_iteration/loop_states_ras.nii.gz"
PYTHONPATH="$ROOT/src:$ROOT/tests" "$ROOT/.venv/bin/python" \
  "$ROOT/tests/java_reference/generate_lesion_repeated_iteration_volume_fixture.py"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
javac -d "$TMP" "$HARNESS"
java -cp "$TMP" LesionToadsRepeatedIterationReferenceHarness "$LOOP_STATES" > "$OUT"
python3 -c 'import hashlib,json,sys; h,o=sys.argv[1:]; r=json.load(open(o)); r["harness_sha256"]=hashlib.sha256(open(h,"rb").read()).hexdigest(); open(o,"w").write(json.dumps(r,indent=2)+"\n")' "$HARNESS" "$OUT"
