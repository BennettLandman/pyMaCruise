#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/classes"
inputs="$repo_root/tests/data/synthetic/folded_cortex"
outputs="$inputs/java_reference"

mkdir -p "$classes" "$outputs"
javac -d "$classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/TGDMReferenceHarness.java"
java -Xmx2g -cp "$classes" TGDMReferenceHarness "$outputs" "$inputs"
(cd "$outputs" && shasum -a 256 *.nii.gz > SHA256SUMS)
