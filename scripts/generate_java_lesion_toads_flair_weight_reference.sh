#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
classes="$repo_root/build/java-reference/lesion-toads-flair-weight-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_flair_weight/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$upstream" "$reference_root/LesionToadsFlairWeightReferenceHarness.java"
java -cp "$classes" LesionToadsFlairWeightReferenceHarness > "$output"
