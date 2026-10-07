#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/tgdm-topology-classes"
central="$repo_root/tests/data/synthetic/tgdm_core_brain_resolution/java_reference/central.nii.gz"
mkdir -p "$classes"
javac -d "$classes" \
  -sourcepath "$reference_root/stubs:$reference_root/upstream" \
  "$reference_root/TGDMCoreReferenceHarness.java" \
  "$reference_root/TGDMTopologyDecisionHarness.java"
java -Xmx2g -cp "$classes" TGDMTopologyDecisionHarness "$central"
