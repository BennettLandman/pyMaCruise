#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-inhomogeneity-degree1-fit"
OUT="$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/java_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsInhomogeneityDegreeOneFitReferenceHarness.java"
java -cp "$BUILD" ToadsInhomogeneityDegreeOneFitReferenceHarness \
  "$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/image.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/memberships.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/best.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/segmentation.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_inhomogeneity_fit_degree1/centroids.txt" > "$OUT"
