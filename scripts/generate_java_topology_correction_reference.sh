#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
case_dir="$repo_root/tests/java_reference/topology_correction"
classes="$repo_root/build/java-reference/topology-correction-classes"
output="$repo_root/build/java-reference/topology-correction-output.txt"
trace="$repo_root/build/java-reference/topology-correction-trace.txt"
mkdir -p "$classes" "$(dirname -- "$output")"

javac -d "$classes" \
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
  "$case_dir/TopologyCorrectionReferenceHarness.java"

if [ "${TOPOLOGY_TRACE:-0}" = 1 ]; then
  fixture=${TOPOLOGY_TRACE_CASE:-spatial_hash}
  java -Dmacruise.topology.trace=true -cp "$classes" TopologyCorrectionReferenceHarness "$fixture" > "$output" 2> "$trace"
  "$repo_root/.venv/bin/python" "$case_dir/convert_output.py" "$fixture" "$output"
else
  for fixture in spatial_hash plateau boundary narrow_branch disconnected_seeds; do
    java -cp "$classes" TopologyCorrectionReferenceHarness "$fixture" > "$output"
    "$repo_root/.venv/bin/python" "$case_dir/convert_output.py" "$fixture" "$output"
  done
fi
