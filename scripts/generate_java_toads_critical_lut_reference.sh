#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-critical-lut-classes"
output="$repo_root/tests/data/synthetic/toads_critical_lut/java_reference.json"
source="${reference_root}/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$source" \
  "$reference_root/ToadsCriticalPointLutReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" \
  "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsCriticalPointLutReferenceHarness > "$output"
