#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-thinning-transition-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_thinning_transition/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java" \
  "$repo_root/tests/java_reference/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$repo_root/tests/java_reference/LesionToadsThinningTransitionReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" LesionToadsThinningTransitionReferenceHarness > "$output"
