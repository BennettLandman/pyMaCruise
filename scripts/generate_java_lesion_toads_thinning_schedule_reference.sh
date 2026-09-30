#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/lesion-thinning-schedule-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_thinning_schedule/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$reference_root/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java" \
  "$reference_root/LesionToadsThinningScheduleReferenceHarness.java"
java -cp "$classes" LesionToadsThinningScheduleReferenceHarness > "$output"
