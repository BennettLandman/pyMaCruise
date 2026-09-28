#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-optimized-similarity-classes"
output="$repo_root/tests/data/synthetic/toads_optimized_similarity/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsOptimizedSimilarityReferenceHarness.java"
java -cp "$classes" ToadsOptimizedSimilarityReferenceHarness > "$output"
