#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workflow_root=${1:-"$repo_root/build/plus-pyramid-source-current/natural-run/plus-surface"}
reference_dir=${2:-"$repo_root/build/plus-tgdm-current-java-2026-10-08"}
report=${3:-"$repo_root/docs/reports/plus-boundary-current-run-java-comparison-2026-10-08.json"}
mkdir -p "$reference_dir/classes"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_tgdm_java_inputs.py" --workflow-root "$workflow_root" --output-dir "$reference_dir"
source_root="$repo_root/tests/java_reference/upstream"
stub_root="$repo_root/tests/java_reference/stubs"
javac -sourcepath "$stub_root:$source_root" -d "$reference_dir/classes" \
  "$source_root/edu/vanderbilt/masi/plugins/CRUISE/utilities/MinFilter3D.java" \
  "$repo_root/tests/java_reference/plus_boundary/PlusBoundarySourceHarness.java"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms128m -Xmx2g -cp "$reference_dir/classes" \
  edu.vanderbilt.masi.plugins.CRUISE.utilities.PlusBoundarySourceHarness \
  "$reference_dir/boundary-labels.bin" "$reference_dir/java-raw-boundary.bin" \
  > "$reference_dir/boundary-java-watchdog.log" 2>&1
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms128m -Xmx2g -cp "$reference_dir/classes" \
  edu.vanderbilt.masi.plugins.CRUISE.utilities.PlusBoundarySourceHarness \
  "$reference_dir/boundary-folded-labels.bin" "$reference_dir/java-folded-raw-boundary.bin" \
  > "$reference_dir/boundary-folded-java-watchdog.log" 2>&1
cat "$reference_dir/boundary-java-watchdog.log"
cat "$reference_dir/boundary-folded-java-watchdog.log"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/compare_plus_boundary_java.py" \
  --reference-dir "$reference_dir" --report "$report"
