#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
input_root="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_inputs"
output_root="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_reference"
classes="$repo_root/build/java-reference/tgdm-core-classes"

mkdir -p "$classes" "$output_root"
javac -d "$classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/TGDMCoreReferenceHarness.java"
java -Xms512m -Xmx8g -cp "$classes" TGDMCoreReferenceHarness "$input_root" "$output_root"
"$repo_root/.venv/bin/python" "$repo_root/scripts/convert_tgdm_core_java_reference.py" \
  --input-root "$input_root" --reference-root "$output_root"
(cd "$output_root" && shasum -a 256 *.nii.gz > SHA256SUMS)
