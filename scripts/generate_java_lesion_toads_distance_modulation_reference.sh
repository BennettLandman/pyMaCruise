#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/lesion-toads-distance-modulation-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_distance_modulation/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$reference_root/LesionToadsDistanceModulationReferenceHarness.java"
java -cp "$classes" LesionToadsDistanceModulationReferenceHarness > "$output"
