#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
NUMERICS="$ROOT/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
BUILD="$ROOT/build/java-reference/toads-default-approx-initial-membership"
OUT="$ROOT/tests/data/synthetic/toads_default_approximate_registration/initial_membership_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" "$NUMERICS" \
  "$ROOT/tests/java_reference/ToadsDefaultApproximateInitialMembershipReferenceHarness.java"
java -cp "$BUILD" ToadsDefaultApproximateInitialMembershipReferenceHarness \
  "$ROOT/tests/data/synthetic/toads_default_approximate_registration/images.nii.gz" \
  "$ROOT/data/atlases/cruise-atlas-10obj-toads2008.txt" > "$OUT"
