#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
fixture_dir="$repo_root/tests/data/synthetic/lesion_data_energy_maps_special"
build_dir="$(mktemp -d "$repo_root/.java-energy-maps-special-build.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT
mkdir -p "$fixture_dir"
javac -d "$build_dir" \
  "$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java" \
  "$repo_root/tests/java_reference/LesionDataEnergyMapsSpecialReferenceHarness.java"
java -cp "$build_dir" LesionDataEnergyMapsSpecialReferenceHarness > "$fixture_dir/java_reference.json"
