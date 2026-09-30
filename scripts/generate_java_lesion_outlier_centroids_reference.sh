#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-outlier-centroids-classes"
output="$repo_root/tests/data/synthetic/lesion_outlier_centroids/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionOutlierCentroidsReferenceHarness.java"
java -cp "$classes" LesionOutlierCentroidsReferenceHarness > "$output"
