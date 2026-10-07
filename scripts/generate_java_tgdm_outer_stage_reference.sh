#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
input_root="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_inputs"
core_root="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_reference"
output_root="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_outer_isolation"
classes="$repo_root/build/java-reference/tgdm-outer-classes"

mkdir -p "$classes" "$output_root"
javac -d "$classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/TGDMCoreReferenceHarness.java" \
  "$reference_root/TGDMOuterStageReferenceHarness.java"
java -Xms512m -Xmx4g -cp "$classes" TGDMOuterStageReferenceHarness \
  "$input_root" "$core_root" "$output_root"
"$repo_root/.venv/bin/python" "$repo_root/scripts/convert_tgdm_core_java_reference.py" \
  --input-root "$input_root" --reference-root "$output_root"
(cd "$output_root" && shasum -a 256 *.nii.gz *.tsv > SHA256SUMS)
