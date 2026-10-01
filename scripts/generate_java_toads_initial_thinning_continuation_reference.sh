#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
fixture="$repo_root/tests/data/synthetic/toads_initial_thinning_startup"
classes="$repo_root/build/java-reference/toads-initial-thinning-continuation-classes"
output="$fixture/continuation_java_reference.json"

mkdir -p "$classes" "$fixture"
javac -d "$classes" \
  "$reference_root/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java" \
  "$reference_root/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java" \
  "$reference_root/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsRelationDiffusionReferenceHarness.java" \
  "$reference_root/ToadsBoundaryDistancesReferenceHarness.java" \
  "$reference_root/ToadsClosestLabelsReferenceHarness.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsInitialThinningStartupReferenceHarness.java" \
  "$reference_root/ToadsInitialStatisticsMembershipContinuationReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsInitialStatisticsMembershipContinuationReferenceHarness "$fixture" > "$output"
