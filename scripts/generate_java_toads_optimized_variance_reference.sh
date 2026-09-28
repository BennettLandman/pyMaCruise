#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-optimized-variance-classes"
output="$repo_root/tests/data/synthetic/toads_optimized_variance/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsOptimizedVarianceReferenceHarness.java"
java -cp "$classes" ToadsOptimizedVarianceReferenceHarness > "$output"
