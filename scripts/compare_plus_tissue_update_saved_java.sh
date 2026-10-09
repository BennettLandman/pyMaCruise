#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workflow_root=${1:-"$repo_root/build/plus-pyramid-source-current/natural-run/plus-surface"}
reference_dir=${2:-"$repo_root/build/plus-tgdm-current-java-2026-10-08"}
report=${3:-"$repo_root/docs/reports/plus-tissue-current-run-java-comparison-2026-10-08.json"}
mkdir -p "$reference_dir/classes"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_topology_java_inputs.py" --workflow-root "$workflow_root" --output-dir "$reference_dir"
source_root="$repo_root/tests/java_reference/upstream"
stub_root="$repo_root/tests/java_reference/stubs"
javac -sourcepath "$stub_root:$source_root" -d "$reference_dir/classes" \
  "$source_root/edu/vanderbilt/masi/plugins/CRUISE/utilities/MASImorph.java" \
  "$source_root/edu/vanderbilt/masi/plugins/CRUISE/MaCRUISE/CombineMAandTOADS.java" \
  "$repo_root/tests/java_reference/plus_tissue_update/PlusTissueUpdateSourceHarness.java"
/usr/bin/java -Xms128m -Xmx2g -cp "$reference_dir/classes" \
  edu.vanderbilt.masi.plugins.CRUISE.MaCRUISE.PlusTissueUpdateSourceHarness \
  "$reference_dir/tissue-arrays.bin" "$reference_dir" > "$reference_dir/tissue-java.log" 2>&1
cat "$reference_dir/tissue-java.log"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/compare_plus_tissue_update_java.py" --reference-dir "$reference_dir" --report "$report"
