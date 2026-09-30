#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" "$repo_root/tests/java_reference/LesionToadsLesionSweepReferenceHarness.java"
java -cp "$classes_dir" LesionToadsLesionSweepReferenceHarness > "$repo_root/tests/data/synthetic/lesion_toads_lesion_sweep/java_reference.json"
java -cp "$classes_dir" LesionToadsLesionSweepReferenceHarness optimized-middle > "$repo_root/tests/data/synthetic/lesion_toads_lesion_sweep/java_reference_optimized_middle.json"
java -cp "$classes_dir" LesionToadsLesionSweepReferenceHarness optimized-last > "$repo_root/tests/data/synthetic/lesion_toads_lesion_sweep/java_reference_optimized_last.json"
