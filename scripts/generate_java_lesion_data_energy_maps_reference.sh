#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
fixture_dir="$repo_root/tests/data/synthetic/lesion_data_energy_maps"
build_dir="$(mktemp -d "$repo_root/.java-energy-maps-build.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT
mkdir -p "$fixture_dir"
javac -d "$build_dir" "$repo_root/tests/java_reference/LesionDataEnergyMapsReferenceHarness.java"
java -cp "$build_dir" LesionDataEnergyMapsReferenceHarness > "$fixture_dir/java_reference.json"
