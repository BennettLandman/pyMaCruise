#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/dura-update-classes"
output="$repo_root/tests/data/synthetic/dura_dynamic_update/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  -sourcepath "$reference_root/dura_stubs:$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/DynamicPressureUpdateReferenceHarness.java"
java -cp "$classes" DynamicPressureUpdateReferenceHarness > "$output"
