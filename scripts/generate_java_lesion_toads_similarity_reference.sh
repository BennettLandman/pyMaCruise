#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
output="$repo_root/tests/data/synthetic/lesion_toads_similarity/java_reference.json"
mkdir -p "$(dirname -- "$output")"
javac -d "$repo_root/build/java-reference/lesion-toads-similarity-classes" \
  "$reference_root/LesionToadsSimilarityReferenceHarness.java"
java -cp "$repo_root/build/java-reference/lesion-toads-similarity-classes" \
  LesionToadsSimilarityReferenceHarness > "$output"
