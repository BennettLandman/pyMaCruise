#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic"
classes="$repo_root/build/java-reference/lesion-topology-thinning-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_topology_thinning/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$upstream/structures/BinaryHeap4D.java" \
  "$upstream/utilities/Numerics.java" \
  "$reference_root/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/LesionToadsTopologyThinningReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" LesionToadsTopologyThinningReferenceHarness > "$output"
