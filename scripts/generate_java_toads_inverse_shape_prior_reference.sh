#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-inverse-shape-prior-classes"
output="$repo_root/tests/data/synthetic/toads_inverse_shape_prior/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsInverseShapePriorReferenceHarness.java"
java -cp "$classes" ToadsInverseShapePriorReferenceHarness > "$output"
