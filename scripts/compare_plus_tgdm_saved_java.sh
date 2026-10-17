#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workflow_root=${1:-"$repo_root/build/plus-pyramid-source-current/natural-run/plus-surface"}
reference_dir=${2:-"$repo_root/build/plus-tgdm-current-java-2026-10-08"}
report=${3:-"$repo_root/docs/reports/plus-tgdm-current-run-java-comparison-2026-10-08.json"}
context_manifest=${4:-}
mkdir -p "$reference_dir/java" "$reference_dir/classes"
if [ -n "$context_manifest" ]; then
  PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_tgdm_java_inputs.py" \
    --workflow-root "$workflow_root" --output-dir "$reference_dir" \
    --context-manifest "$context_manifest" --allow-python-stage-differences
else
  PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_tgdm_java_inputs.py" \
    --workflow-root "$workflow_root" --output-dir "$reference_dir"
fi
reference_root="$repo_root/tests/java_reference"
javac -d "$reference_dir/classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/PlusTGDMReferenceHarness.java"
"$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g -cp "$reference_dir/classes" \
  PlusTGDMReferenceHarness "$reference_dir/arrays.bin" "$reference_dir/java" \
  precomputed_initial_level_set "$reference_dir/brain-mask.bin" \
  > "$reference_dir/java-watchdog.log" 2>&1
cat "$reference_dir/java-watchdog.log"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/compare_plus_tgdm_java.py" \
  --workflow-root "$workflow_root" --reference-dir "$reference_dir" --report "$report"
