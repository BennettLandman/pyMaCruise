#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT
javac -d "$build_dir" "$repo_root/tests/java_reference/LesionToadsStandard13VolumeReferenceHarness.java"
java -cp "$build_dir" LesionToadsStandard13VolumeReferenceHarness \
  > "$repo_root/tests/data/synthetic/lesion_standard13_volume/java_reference.json"
