#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/classes"
outputs="$repo_root/tests/data/synthetic/tgdm/java_reference"

mkdir -p "$classes" "$outputs"
javac -d "$classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/TGDMReferenceHarness.java"
java -cp "$classes" TGDMReferenceHarness "$outputs"
