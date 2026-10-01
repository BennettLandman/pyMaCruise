#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-default-approximate-centroid"
FIXTURE="$ROOT/tests/data/synthetic/toads_approximate_centroids"
OUT="$FIXTURE/java_reference.json"
POWER="$ROOT/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
mkdir -p "$BUILD" "$FIXTURE"
javac -d "$BUILD" "$POWER" "$ROOT/tests/java_reference/ToadsDefaultApproximateCentroidReferenceHarness.java"
java -cp "$BUILD" ToadsDefaultApproximateCentroidReferenceHarness \
  "$FIXTURE/images.nii.gz" \
  "$FIXTURE/field.nii.gz" \
  "$FIXTURE/memberships.nii.gz" \
  "$FIXTURE/best.nii.gz" \
  "$FIXTURE/centroids.txt" \
  "$FIXTURE/optimized_exponents.txt" > "$OUT"
