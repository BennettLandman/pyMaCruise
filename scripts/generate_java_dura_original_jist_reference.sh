#!/bin/sh
set -eu

if [ "${JIST_TEST_JAR:-}" = "" ] || [ "${JIST_LIB_DIR:-}" = "" ]; then
  echo "Set JIST_TEST_JAR and JIST_LIB_DIR to the supplied JIST build." >&2
  exit 2
fi

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/dura-full-evolution-classes"
runtime="$repo_root/build/java-reference/dura-runtime-original-jist"
output=${DURA_TGDM_OUTPUT:-"$repo_root/tests/data/synthetic/dura_tgdm_full_evolution/java_reference.json"}
nx=${DURA_TGDM_NX:-11}
ny=${DURA_TGDM_NY:-10}
nz=${DURA_TGDM_NZ:-9}
folded=${DURA_TGDM_FOLDED:-false}
case "$output" in
  /*) ;;
  *) output="$repo_root/$output" ;;
esac

"$repo_root/scripts/generate_java_dura_tgdm_full_evolution_reference.sh"
python3 "$repo_root/scripts/prepare_java_dura_original_jist_runtime.py" \
  "$JIST_TEST_JAR" "$classes" "$runtime" >/dev/null
mkdir -p "$(dirname -- "$output")"
java -Ddura.original_jist=true -Ddura.nx="$nx" -Ddura.ny="$ny" \
  -Ddura.nz="$nz" -Ddura.folded="$folded" \
  -cp "$runtime:$runtime/../jist-runtime-without-selected-conflicts.jar:$JIST_LIB_DIR/vecmath-1.3.1.jar:$JIST_LIB_DIR/commons-math3-3.0.jar" \
  DuraTGDMFullEvolutionReferenceHarness "$output" >/dev/null
