#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
classes_dir="$repo_root/build/java-reference/toads-object-mask-outlier-sweep"
mkdir -p "$classes_dir"
javac -d "$classes_dir" "$repo_root/tests/java_reference/ToadsObjectMaskOutlierSweepReferenceHarness.java"
java -cp "$classes_dir" ToadsObjectMaskOutlierSweepReferenceHarness > "$repo_root/tests/data/synthetic/toads_object_mask_outlier_sweep/java_reference.json"
