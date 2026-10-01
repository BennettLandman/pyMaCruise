#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-default-exact-centroids-classes"
output="$repo_root/tests/data/synthetic/toads_default_first_thinning/centroid_reference.json"
power_source="$reference_root/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" \
  "$reference_root/ToadsDefaultExactCentroidReferenceHarness.java"
java -cp "$classes" ToadsDefaultExactCentroidReferenceHarness > "$output"
