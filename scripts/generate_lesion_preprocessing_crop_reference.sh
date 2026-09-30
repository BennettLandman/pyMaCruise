#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
reference_dir="$(mktemp -d)"
trap 'rm -rf "$reference_dir"' EXIT
javac -d "$reference_dir" "$repo_root/tests/java_reference/LesionPreprocessingCropReferenceHarness.java"
java -cp "$reference_dir" LesionPreprocessingCropReferenceHarness > "$repo_root/tests/data/synthetic/lesion_preprocessing_crop/java_reference.json"
