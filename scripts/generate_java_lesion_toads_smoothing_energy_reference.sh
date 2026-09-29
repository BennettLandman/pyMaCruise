#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/lesion-toads-smoothing-energy-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_smoothing_energy/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$reference_root/LesionToadsSmoothingEnergyReferenceHarness.java"
java -cp "$classes" LesionToadsSmoothingEnergyReferenceHarness > "$output"
