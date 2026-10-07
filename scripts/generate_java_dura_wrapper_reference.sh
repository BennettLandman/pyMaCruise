#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/dura-attached-no-dura-classes"
output=${DURA_WRAPPER_OUTPUT:-"$repo_root/build/java-reference/dura-attached-no-dura/reference.json"}
nx=${DURA_WRAPPER_NX:-17}
ny=${DURA_WRAPPER_NY:-17}
nz=${DURA_WRAPPER_NZ:-17}
heap=${DURA_WRAPPER_JAVA_HEAP:-2g}

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  -sourcepath "$reference_root/topology_correction/upstream:$reference_root/topology_correction/stubs:$reference_root/upstream:$reference_root/stubs:$reference_root/dura_stubs" \
  "$reference_root/DuraAttachedNoDuraReferenceHarness.java"
java -Xmx"$heap" -Ddura.nx="$nx" -Ddura.ny="$ny" -Ddura.nz="$nz" \
  -cp "$classes" DuraAttachedNoDuraReferenceHarness "$output"
