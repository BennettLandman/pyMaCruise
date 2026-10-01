#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-default-exact-variances-classes"
output="$repo_root/tests/data/synthetic/toads_default_first_thinning/variance_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$reference_root/ToadsDefaultExactVarianceReferenceHarness.java"
java -cp "$classes" ToadsDefaultExactVarianceReferenceHarness > "$output"
