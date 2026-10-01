#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-default-11-membership-continuation-classes"
output="$repo_root/tests/data/synthetic/toads_default_11_continue_thinning/membership_reference.json"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/structures"
power_source="$reference_root/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" \
  "$upstream/CriticalPointLUT.java" \
  "$upstream/BinaryHeap4D.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsDefault11ContinueThinningReferenceHarness.java" \
  "$reference_root/ToadsDefault11MembershipContinuationReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsDefault11MembershipContinuationReferenceHarness > "$output"
python3 - "$reference_root/ToadsDefault11MembershipContinuationReferenceHarness.java" "$output" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

harness = Path(sys.argv[1])
reference = Path(sys.argv[2])
payload = json.loads(reference.read_text(encoding="utf-8"))
payload["harness_sha256"] = hashlib.sha256(harness.read_bytes()).hexdigest()
reference.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
PY
