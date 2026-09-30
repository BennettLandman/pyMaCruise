#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
heap_source="$reference_root/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"
classes="$repo_root/build/java-reference/lesion-thinning-queue-run-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_thinning_queue_run/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$heap_source" "$reference_root/LesionToadsThinningQueueRunReferenceHarness.java"
java -cp "$classes" LesionToadsThinningQueueRunReferenceHarness > "$output"
