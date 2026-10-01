#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-default-approximate-membership-pass"
FIXTURE="$ROOT/tests/data/synthetic/toads_default_approximate_pass"
mkdir -p "$BUILD" "$FIXTURE"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsDefaultApproximateMembershipPassHarness.java"
java -cp "$BUILD" ToadsDefaultApproximateMembershipPassHarness \
  "$FIXTURE/images.nii.gz" "$FIXTURE/field.nii.gz" \
  "$FIXTURE/memberships.nii.gz" "$FIXTURE/best.nii.gz" \
  "$FIXTURE/transformed_shapes.nii.gz" "$FIXTURE/direct_relation.nii.gz" \
  "$FIXTURE/centroids.txt" "$FIXTURE/standard_deviations.txt" \
  "$FIXTURE/variance_weights.txt" "$FIXTURE/exponents.txt" \
  "$FIXTURE/modality.txt" "$FIXTURE/smoothing_priors.txt" \
  "$FIXTURE/shape_priors.txt" "$FIXTURE/parameters.txt" \
  > "$FIXTURE/java_reference.json"
