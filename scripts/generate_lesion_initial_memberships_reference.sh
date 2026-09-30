#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
reference_dir="$(mktemp -d)"
trap 'rm -rf "$reference_dir"' EXIT
javac -d "$reference_dir" "$repo_root/tests/java_reference/LesionInitialMembershipsReferenceHarness.java"
java -cp "$reference_dir" LesionInitialMembershipsReferenceHarness > "$repo_root/tests/data/synthetic/lesion_initial_memberships/java_reference.json"
