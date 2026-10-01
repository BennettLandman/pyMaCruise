#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-default-approximate-variance"
FIXTURE="$ROOT/tests/data/synthetic/toads_approximate_variances"
OUT="$FIXTURE/java_reference.json"
mkdir -p "$BUILD" "$FIXTURE"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsDefaultApproximateVarianceReferenceHarness.java"
java -cp "$BUILD" ToadsDefaultApproximateVarianceReferenceHarness \
  "$FIXTURE/images.nii.gz" \
  "$FIXTURE/field.nii.gz" \
  "$FIXTURE/memberships.nii.gz" \
  "$FIXTURE/best.nii.gz" \
  "$FIXTURE/segmentation.nii.gz" \
  "$FIXTURE/updated_centroids.txt" > "$OUT"
