#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
source_file="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
expected_sha="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
actual_sha=$(shasum -a 256 "$source_file" | awk '{print $1}')
if [ "$actual_sha" != "$expected_sha" ]; then
    echo "Unexpected MultipleToadSegmentation source hash: $actual_sha" >&2
    exit 1
fi
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" \
    "$repo_root/tests/java_reference/ToadsRelationListReferenceHarness.java" \
    "$repo_root/tests/java_reference/ToadsRelationDiffusionReferenceHarness.java" \
    "$repo_root/tests/java_reference/ToadsBoundaryDistancesReferenceHarness.java" \
    "$repo_root/tests/java_reference/ToadsClosestLabelsReferenceHarness.java" \
    "$repo_root/tests/java_reference/ToadsInitialRelationStateReferenceHarness.java"
java -cp "$classes_dir" ToadsInitialRelationStateReferenceHarness > \
    "$repo_root/tests/data/synthetic/toads_initial_relation_state/java_reference.json"
