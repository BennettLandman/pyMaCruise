#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-shape-relations-classes"
output="$repo_root/tests/data/synthetic/toads_shape_relations/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsShapeRelationReferenceHarness.java"
java -cp "$classes" ToadsShapeRelationReferenceHarness > "$output"
