#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-weighted-centroid-classes"
output="$repo_root/tests/data/synthetic/lesion_constrained_centroids/weighted_java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionWeightedCentroidsReferenceHarness.java"
java -cp "$classes" LesionWeightedCentroidsReferenceHarness > "$output"
