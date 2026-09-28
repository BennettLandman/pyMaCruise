#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-shape-bounds-classes"
output="$repo_root/tests/data/synthetic/toads_shape_bounds/java_reference.json"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$reference_root/stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$upstream/utilities/Numerics.java" \
  "$upstream/structures/RotationMatrix.java" \
  "$upstream/methods/ParametricTransform.java" \
  "$upstream/libraries/ImageFunctions.java" \
  "$reference_root/ToadsShapeBoundsReferenceHarness.java"
java -cp "$classes" ToadsShapeBoundsReferenceHarness > "$output"
