#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/structures"
classes="$repo_root/build/java-reference/lesion-composed-thinning-event-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_composed_thinning_event/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/BinaryHeap4D.java" \
  "$upstream/CriticalPointLUT.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/LesionToadsComposedThinningEventReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" LesionToadsComposedThinningEventReferenceHarness > "$output"
