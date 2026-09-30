#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$repo_root/build/java-reference/toads-pair-similarity"
mkdir -p "$build_dir"
mkdir -p "$repo_root/tests/data/synthetic/toads_pair_similarity"
javac -d "$build_dir" "$repo_root/tests/java_reference/ToadsPairSimilarityReferenceHarness.java"
java -cp "$build_dir" ToadsPairSimilarityReferenceHarness \
  > "$repo_root/tests/data/synthetic/toads_pair_similarity/java_reference.json"
