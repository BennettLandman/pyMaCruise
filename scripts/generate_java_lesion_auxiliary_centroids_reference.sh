#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-auxiliary-centroids-classes"
output="$repo_root/tests/data/synthetic/lesion_special_centroids/java_reference.json"
power_source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" \
  "$repo_root/tests/java_reference/LesionAuxiliaryCentroidsReferenceHarness.java"
java -cp "$classes" LesionAuxiliaryCentroidsReferenceHarness > "$output"
