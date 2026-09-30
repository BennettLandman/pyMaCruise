#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-basic-centroids-classes"
output="$repo_root/tests/data/synthetic/lesion_variances/basic_centroids_java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionBasicCentroidsReferenceHarness.java"
java -cp "$classes" LesionBasicCentroidsReferenceHarness > "$output"
