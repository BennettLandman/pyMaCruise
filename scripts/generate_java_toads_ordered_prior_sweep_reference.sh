#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-ordered-prior-sweep-classes"
output="$repo_root/tests/data/synthetic/toads_ordered_prior_sweep/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsOrderedPriorSweepReferenceHarness.java"
java -cp "$classes" ToadsOrderedPriorSweepReferenceHarness > "$output"
