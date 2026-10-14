#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/dura-base-input-classes"
input=${1:-"$repo_root/build/java-reference/dura-base-input.bin"}
output=${2:-"$repo_root/build/java-reference/dura-base-output"}
heap=${DURA_WRAPPER_JAVA_HEAP:-8g}
nx=${DURA_WRAPPER_NX:-172}
ny=${DURA_WRAPPER_NY:-216}
nz=${DURA_WRAPPER_NZ:-164}
iterations=${DURA_WRAPPER_ITERATIONS:-1}
mkdir -p "$classes" "$output"
javac -d "$classes" \
  -sourcepath "$reference_root/topology_correction/upstream:$reference_root/topology_correction/stubs:$reference_root/upstream:$reference_root/stubs:$reference_root/dura_stubs" \
  "$reference_root/DuraAttachedNoDuraReferenceHarness.java" \
  "$reference_root/DuraBaseInputsReferenceHarness.java"
java -Xmx"$heap" -Ddura.nx="$nx" -Ddura.ny="$ny" -Ddura.nz="$nz" \
  -Ddura.iterations="$iterations" -Ddura.membership="${DURA_WRAPPER_MEMBERSHIP:-0.4}" -cp "$classes" \
  DuraBaseInputsReferenceHarness "$input" "$output"
