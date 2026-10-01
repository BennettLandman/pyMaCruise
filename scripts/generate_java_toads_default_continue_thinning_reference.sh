#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-default-continue-thinning-classes"
output="$repo_root/tests/data/synthetic/toads_default_continue_thinning/java_reference.json"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/structures"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/CriticalPointLUT.java" \
  "$upstream/BinaryHeap4D.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsRelationDiffusionReferenceHarness.java" \
  "$reference_root/ToadsBoundaryDistancesReferenceHarness.java" \
  "$reference_root/ToadsClosestLabelsReferenceHarness.java" \
  "$reference_root/ToadsGrowingPassReferenceHarness.java" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsDefaultContinueThinningReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsDefaultContinueThinningReferenceHarness > "$output"
