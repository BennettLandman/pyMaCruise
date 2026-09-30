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
