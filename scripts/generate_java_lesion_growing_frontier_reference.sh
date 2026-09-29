#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
classes="$repo_root/build/java-reference/lesion-growing-frontier-classes"
output="$repo_root/tests/data/synthetic/lesion_growing_frontier/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/utilities/Numerics.java" \
  "$upstream/structures/BinaryHeap4D.java" \
  "$reference_root/LesionGrowingFrontierReferenceHarness.java"
java -cp "$classes" LesionGrowingFrontierReferenceHarness > "$output"
