#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
HARNESS="$ROOT/tests/java_reference/LesionRelationAwareMembershipReferenceHarness.java"
OUT="$ROOT/tests/data/synthetic/lesion_continuing_iteration/relation_dispatch_java_reference.json"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
javac -d "$TMP" "$HARNESS"
java -cp "$TMP" LesionRelationAwareMembershipReferenceHarness > "$OUT"
python3 -c 'import hashlib,json,sys; h,o=sys.argv[1:]; r=json.load(open(o)); r["harness_sha256"]=hashlib.sha256(open(h,"rb").read()).hexdigest(); open(o,"w").write(json.dumps(r,separators=(",",":"))+"\n")' "$HARNESS" "$OUT"
