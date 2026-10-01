#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-inhomogeneity-transition"
FIXTURE="$ROOT/tests/data/synthetic/toads_inhomogeneity_transition"
OUT="$FIXTURE/java_reference.json"
mkdir -p "$BUILD" "$FIXTURE"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsInhomogeneityTransitionReferenceHarness.java"
java -cp "$BUILD" ToadsInhomogeneityTransitionReferenceHarness \
  "$FIXTURE/images.nii.gz" \
  "$FIXTURE/memberships.nii.gz" \
  "$FIXTURE/best.nii.gz" \
  "$FIXTURE/segmentation.nii.gz" \
  "$FIXTURE/attached_fields.nii.gz" \
  "$FIXTURE/centroids.txt" > "$OUT"
