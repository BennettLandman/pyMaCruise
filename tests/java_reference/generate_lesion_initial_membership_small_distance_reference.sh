#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
FIXTURE="$ROOT/tests/data/synthetic/lesion_initial_membership_small_distance"
HARNESS="$ROOT/tests/java_reference/LesionInitialMembershipSmallDistanceReferenceHarness.java"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
"$ROOT/.venv/bin/python" "$FIXTURE/generate_fixture.py"
javac -d "$TMP" "$HARNESS"
java -cp "$TMP" LesionInitialMembershipSmallDistanceReferenceHarness > "$FIXTURE/java_reference.json"
python3 -c 'import hashlib,json,sys; h,o=sys.argv[1:]; r=json.load(open(o)); r["harness_sha256"]=hashlib.sha256(open(h,"rb").read()).hexdigest(); open(o,"w").write(json.dumps(r,separators=(",",":"))+"\n")' "$HARNESS" "$FIXTURE/java_reference.json"
