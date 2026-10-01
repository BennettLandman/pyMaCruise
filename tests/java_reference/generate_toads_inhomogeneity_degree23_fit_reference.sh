#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-inhomogeneity-degree23-fit"
FIXTURE="$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree23"
OUT="$FIXTURE/java_reference.json"
mkdir -p "$BUILD" "$FIXTURE"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsInhomogeneityDegree23FitReferenceHarness.java"
java -cp "$BUILD" ToadsInhomogeneityDegree23FitReferenceHarness \
  "$FIXTURE/image.nii.gz" \
  "$FIXTURE/memberships.nii.gz" \
  "$FIXTURE/best.nii.gz" \
  "$FIXTURE/segmentation.nii.gz" \
  "$FIXTURE/centroids.txt" > "$OUT"
