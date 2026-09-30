#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
classes="$repo_root/build/java-reference/lesion-thinning-pass-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_thinning_pass/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/structures/BinaryHeap4D.java" \
  "$upstream/utilities/Numerics.java" \
  "$reference_root/LesionToadsThinningPassReferenceHarness.java"
java -cp "$classes" LesionToadsThinningPassReferenceHarness > "$output"
