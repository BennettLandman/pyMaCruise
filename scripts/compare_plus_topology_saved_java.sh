#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workflow_root=${1:-"$repo_root/build/plus-pyramid-source-current/natural-run/plus-surface"}
reference_dir=${2:-"$repo_root/build/plus-tgdm-current-java-2026-10-08"}
report=${3:-"$repo_root/docs/reports/plus-topology-current-run-java-comparison-2026-10-08.json"}
mkdir -p "$reference_dir/classes"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/prepare_plus_topology_java_inputs.py" --workflow-root "$workflow_root" --output-dir "$reference_dir"
case_dir="$repo_root/tests/java_reference/topology_correction"
javac -d "$reference_dir/classes" \
  "$case_dir/stubs/gov/nih/mipav/model/algorithms/AlgorithmBase.java" \
  "$case_dir/stubs/gov/nih/mipav/model/structures/ModelImage.java" \
  "$case_dir/stubs/gov/nih/mipav/view/ViewUserInterface.java" \
  "$case_dir/stubs/gov/nih/mipav/view/ViewJProgressBar.java" \
  "$case_dir/stubs/gov/nih/mipav/view/Preferences.java" \
  "$case_dir/stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$case_dir/stubs/edu/jhmi/rad/medic/libraries/ObjectProcessing.java" \
  "$case_dir/stubs/WildMagic/LibFoundation/Mathematics/Vector3f.java" \
  "$case_dir/upstream/edu/jhmi/rad/medic/structures/BinaryTree.java" \
  "$case_dir/upstream/edu/jhmi/rad/medic/libraries/NeighborConnectivity.java" \
  "$case_dir/upstream/edu/jhmi/rad/medic/libraries/TopologyPropagation.java" \
  "$case_dir/upstream/edu/jhmi/rad/medic/algorithms/AlgorithmTopologyCorrection.java" \
  "$case_dir/PlusTopologyReferenceHarness.java"
"$repo_root/.venv/bin/python" "$repo_root/scripts/run_with_rss_watchdog.py" --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms128m -Xmx6g -cp "$reference_dir/classes" PlusTopologyReferenceHarness \
  "$reference_dir/topology-arrays.bin" "$reference_dir/topology-java.raw.gz" > "$reference_dir/topology-java-watchdog.log" 2>&1
cat "$reference_dir/topology-java-watchdog.log"
PYTHONPATH="$repo_root/src" "$repo_root/.venv/bin/python" "$repo_root/scripts/compare_plus_topology_java.py" --reference-dir "$reference_dir" --report "$report"
