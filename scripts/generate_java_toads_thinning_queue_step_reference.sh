#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-thinning-queue-step-classes"
output="$repo_root/tests/data/synthetic/toads_thinning_queue_step/java_reference.json"
lut_source="$reference_root/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java"
heap_source="$reference_root/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$lut_source" "$heap_source" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsThinningQueueStepReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsThinningQueueStepReferenceHarness > "$output"
