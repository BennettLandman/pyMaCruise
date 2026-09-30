#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
vandy_root="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src"
toads_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
power_source="$vandy_root/edu/jhmi/rad/medic/utilities/PowerTable.java"
expected_toads="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
expected_power="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
actual_toads=$(shasum -a 256 "$toads_source" | awk '{print $1}')
actual_power=$(shasum -a 256 "$power_source" | awk '{print $1}')
if [ "$actual_toads" != "$expected_toads" ]; then
    echo "Unexpected MultipleToadSegmentation source hash: $actual_toads" >&2
    exit 1
fi
if [ "$actual_power" != "$expected_power" ]; then
    echo "Unexpected PowerTable source hash: $actual_power" >&2
    exit 1
fi
classes_dir="$repo_root/build/java-reference/toads-optimized-dispatch"
mkdir -p "$classes_dir"
javac -d "$classes_dir" "$power_source" "$repo_root/tests/java_reference/ToadsOptimizedMixedDispatchReferenceHarness.java"
java -cp "$classes_dir" ToadsOptimizedMixedDispatchReferenceHarness \
    "$repo_root/tests/data/synthetic/toads_optimized_dispatch/images.nii.gz" \
    "$repo_root/tests/data/synthetic/toads_optimized_dispatch/field.nii.gz" \
    "$repo_root/tests/data/synthetic/toads_optimized_dispatch/psf_bias.nii.gz" \
    "$repo_root/tests/data/synthetic/toads_optimized_dispatch/folded_ribbon.nii.gz" \
    > "$repo_root/tests/data/synthetic/toads_optimized_dispatch/java_reference.json"
