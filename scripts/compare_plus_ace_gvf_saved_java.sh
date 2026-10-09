#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workflow_root=${1:-"$repo_root/build/plus-pyramid-source-current/natural-run/plus-surface"}
reference_dir=${2:-"$repo_root/build/plus-tgdm-current-java-2026-10-08"}
report=${3:-"$repo_root/docs/reports/plus-ace-gvf-current-run-java-comparison-2026-10-08.json"}
mkdir -p "$reference_dir/classes" "$reference_dir/ace-intermediate"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_tgdm_java_inputs.py" --workflow-root "$workflow_root" --output-dir "$reference_dir"
source_root="$repo_root/tests/java_reference"
javac -d "$reference_dir/classes" -sourcepath "$source_root/stubs:$source_root/upstream" "$source_root/ACEReferenceHarness.java" "$source_root/GVFReferenceHarness.java"
"$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms128m -Xmx2g -cp "$reference_dir/classes" ACEReferenceHarness \
  "$reference_dir/ace-input.bin" "$reference_dir/ace-output.bin" "$reference_dir/ace-intermediate" 3 > "$reference_dir/ace-java-watchdog.log" 2>&1
cat "$reference_dir/ace-java-watchdog.log"
"$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms128m -Xmx2g -cp "$reference_dir/classes" GVFReferenceHarness \
  "$reference_dir/gvf-input.bin" "$reference_dir/gvf-output.bin" 0.2 200 > "$reference_dir/gvf-java-watchdog.log" 2>&1
cat "$reference_dir/gvf-java-watchdog.log"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/compare_plus_ace_gvf_java.py" --reference-dir "$reference_dir" --report "$report"
