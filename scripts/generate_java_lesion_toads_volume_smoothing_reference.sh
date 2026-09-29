#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" "$repo_root/tests/java_reference/LesionToadsVolumeSmoothingReferenceHarness.java"
java -cp "$classes_dir" LesionToadsVolumeSmoothingReferenceHarness > "$repo_root/tests/data/synthetic/lesion_toads_volume_smoothing/java_reference.json"
