#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-mixed-centroids-classes"
output="$repo_root/tests/data/synthetic/lesion_mixed_centroids/java_reference.json"
power_source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" \
  "$repo_root/tests/java_reference/LesionMixedCentroidsReferenceHarness.java"
java -cp "$classes" LesionMixedCentroidsReferenceHarness > "$output"
