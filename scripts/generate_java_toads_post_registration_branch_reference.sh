#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
source_file="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java"
expected_sha="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
actual_sha=$(shasum -a 256 "$source_file" | awk '{print $1}')
if [ "$actual_sha" != "$expected_sha" ]; then
    echo "Unexpected AlgorithmMultipleToads source hash: $actual_sha" >&2
    exit 1
fi
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" "$repo_root/tests/java_reference/ToadsPostRegistrationBranchReferenceHarness.java"
java -cp "$classes_dir" ToadsPostRegistrationBranchReferenceHarness > "$repo_root/tests/data/synthetic/toads_post_registration_branch/java_reference.json"
