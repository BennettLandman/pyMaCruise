#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-p18czo-nonzero-shapes-next-growth-classes"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
mkdir -p "$classes"
javac -d "$classes" \
  "$upstream/structures/BinaryHeap4D.java" \
  "$upstream/structures/CriticalPointLUT.java" \
  "$upstream/utilities/PowerTable.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsDefault11ContinueThinningReferenceHarness.java" \
  "$reference_root/ToadsP18czoNonzeroShapesNextGrowthReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsP18czoNonzeroShapesNextGrowthReferenceHarness "$@"
