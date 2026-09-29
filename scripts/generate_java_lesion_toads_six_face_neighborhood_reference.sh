#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" "$repo_root/tests/java_reference/LesionToadsSixFaceNeighborhoodReferenceHarness.java"
java -cp "$classes_dir" LesionToadsSixFaceNeighborhoodReferenceHarness > "$repo_root/tests/data/synthetic/lesion_toads_six_face_neighborhood/java_reference.json"
