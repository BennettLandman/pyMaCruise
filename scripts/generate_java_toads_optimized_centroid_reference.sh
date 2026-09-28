#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-optimized-centroid-classes"
output="$repo_root/tests/data/synthetic/toads_optimized_centroid/java_reference.json"
power_source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" "$repo_root/tests/java_reference/ToadsOptimizedCentroidReferenceHarness.java"
java -cp "$classes" ToadsOptimizedCentroidReferenceHarness > "$output"
