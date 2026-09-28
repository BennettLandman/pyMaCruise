#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-shape-weighted-membership-classes"
output="$repo_root/tests/data/synthetic/toads_shape_weighted_membership/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsShapeWeightedMembershipReferenceHarness.java"
java -cp "$classes" ToadsShapeWeightedMembershipReferenceHarness > "$output"
