#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-default-11-continue-thinning-classes"
output="$repo_root/tests/data/synthetic/toads_default_11_continue_thinning/java_reference.json"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/structures"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/CriticalPointLUT.java" \
  "$upstream/BinaryHeap4D.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsDefault11ContinueThinningReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsDefault11ContinueThinningReferenceHarness > "$output"
