#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/dura-full-evolution-classes"
output=${DURA_TGDM_OUTPUT:-"$repo_root/tests/data/synthetic/dura_tgdm_full_evolution/java_reference.json"}
nx=${DURA_TGDM_NX:-11}
ny=${DURA_TGDM_NY:-10}
nz=${DURA_TGDM_NZ:-9}
folded=${DURA_TGDM_FOLDED:-false}
upstream="$reference_root/upstream/edu/jhu/ece/iacl/algorithms"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  -sourcepath "$reference_root/upstream:$reference_root/stubs:$reference_root/dura_stubs" \
  "$upstream/tgdm/GenericTGDM.java" \
  "$upstream/tgdm/DynamicPressureTGDM.java" \
  "$upstream/tgdm/DynamicPressureCortexProfile.java" \
  "$upstream/tgdm/ProfileTGDM.java" \
  "$upstream/gvf/FastMarchingGradient.java" \
  "$reference_root/DuraTGDMFullEvolutionReferenceHarness.java"
java -Ddura.nx="$nx" -Ddura.ny="$ny" -Ddura.nz="$nz" \
  -Ddura.folded="$folded" -cp "$classes" \
  DuraTGDMFullEvolutionReferenceHarness > "$output"
