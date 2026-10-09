#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture="$repo_root/tests/data/synthetic/labeled_folded_cortex"
build_root="$repo_root/build/labeled-folded-cortex-source-selected"
reference_root="$fixture/java_reference"
classes="$build_root/java_classes"
mkdir -p "$reference_root" "$classes"

javac -d "$classes" \
  -sourcepath "$repo_root/tests/java_reference/stubs:$repo_root/tests/java_reference/upstream" \
  "$repo_root/tests/java_reference/TGDMReferenceHarness.java"

"$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" \
  --limit-gb 8 --interval-seconds 0.2 --report-interval-seconds 30 -- \
  java -Xmx4g -cp "$classes" TGDMReferenceHarness "$reference_root" "$fixture" 0.6

(cd "$reference_root" && shasum -a 256 *.nii.gz > SHA256SUMS)
