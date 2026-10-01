#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT
source_root="$repo_root/tests/java_reference/upstream"
javac -d "$tmp_dir" \
  "$source_root/edu/jhmi/rad/medic/utilities/PowerTable.java" \
  "$repo_root/tests/java_reference/ToadsPostThinningExactIterationReferenceHarness.java"
java -cp "$tmp_dir" ToadsPostThinningExactIterationReferenceHarness "$@"
