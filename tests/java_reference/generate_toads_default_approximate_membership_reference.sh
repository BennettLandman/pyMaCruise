#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
UPSTREAM="$ROOT/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
BUILD="$ROOT/build/java-reference/toads-default-approximate-membership"
OUT="$ROOT/tests/data/synthetic/toads_default_approximate_membership/java_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" "$UPSTREAM" \
  "$ROOT/tests/java_reference/ToadsDefaultApproximateMembershipReferenceHarness.java"
java -cp "$BUILD" ToadsDefaultApproximateMembershipReferenceHarness \
  "$ROOT/tests/data/synthetic/toads_default_approximate_membership/images.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_default_approximate_membership/transformed_shapes.nii.gz" \
  "$ROOT/data/atlases/cruise-atlas-10obj-toads2008.txt" > "$OUT"
python3 - "$ROOT/tests/java_reference/ToadsDefaultApproximateMembershipReferenceHarness.java" "$OUT" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

harness = Path(sys.argv[1])
output = Path(sys.argv[2])
reference = json.loads(output.read_text(encoding="utf-8"))
reference["harness_sha256"] = hashlib.sha256(harness.read_bytes()).hexdigest()
output.write_text(json.dumps(reference, separators=(",", ":")) + "\n", encoding="utf-8")
PY
python3 -c 'import hashlib,json,sys; out,*paths=sys.argv[1:]; r=json.load(open(out)); names=("harness","numerics","images","transformed_shapes","descriptor"); r["p18czo_evidence_sha256"]={n:hashlib.sha256(open(p,"rb").read()).hexdigest() for n,p in zip(names,paths)}; open(out,"w").write(json.dumps(r,separators=(",",":"))+"\n")' \
  "$OUT" \
  "$ROOT/tests/java_reference/ToadsDefaultApproximateMembershipReferenceHarness.java" \
  "$UPSTREAM" \
  "$ROOT/tests/data/synthetic/toads_default_approximate_membership/images.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_default_approximate_membership/transformed_shapes.nii.gz" \
  "$ROOT/data/atlases/cruise-atlas-10obj-toads2008.txt"
