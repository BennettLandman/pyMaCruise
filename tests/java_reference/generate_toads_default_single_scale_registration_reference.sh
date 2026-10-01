#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
UPSTREAM="$ROOT/tests/java_reference/upstream/edu/jhmi/rad/medic"
BUILD="$ROOT/build/java-reference/toads-default-single-scale-registration"
OUT="$ROOT/tests/data/synthetic/toads_default_single_scale_registration/java_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" \
  "$UPSTREAM/utilities/MedicUtilPublic.java" \
  "$UPSTREAM/utilities/Numerics.java" \
  "$UPSTREAM/structures/RotationMatrix.java" \
  "$UPSTREAM/methods/ParametricTransform.java" \
  "$UPSTREAM/libraries/ImageFunctions.java" \
  "$ROOT/tests/java_reference/ToadsDefaultSingleScaleRegistrationReferenceHarness.java"
java -cp "$BUILD" ToadsDefaultSingleScaleRegistrationReferenceHarness \
  "$ROOT/tests/data/synthetic/toads_default_single_scale_registration/images.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_default_single_scale_registration/atlas_shapes.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_default_single_scale_registration/transformed_shapes.nii.gz" > "$OUT"
