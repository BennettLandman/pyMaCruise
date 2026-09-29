#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-binary-heap-classes"
output="$repo_root/tests/data/synthetic/toads_binary_heap/java_reference.json"
heap_source="$reference_root/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$heap_source" \
  "$reference_root/ToadsBinaryHeap4DReferenceHarness.java"
java -cp "$classes" ToadsBinaryHeap4DReferenceHarness > "$output"
