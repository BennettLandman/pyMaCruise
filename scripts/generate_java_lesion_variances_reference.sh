#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-variance-classes"
output="$repo_root/tests/data/synthetic/lesion_variances/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionVariancesReferenceHarness.java"
java -cp "$classes" LesionVariancesReferenceHarness > "$output"
