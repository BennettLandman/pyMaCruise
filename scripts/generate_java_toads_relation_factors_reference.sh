#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
source_root="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic"
check_source() {
    file=$1
    expected=$2
    actual=$(shasum -a 256 "$file" | awk '{print $1}')
    if [ "$actual" != "$expected" ]; then
        echo "Unexpected source hash for $file: $actual" >&2
        exit 1
    fi
}
check_source "$source_root/methods/MultipleToadSegmentation.java" \
    "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
check_source "$source_root/utilities/Numerics.java" \
    "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
classes_dir=$(mktemp -d)
trap 'rm -rf "$classes_dir"' EXIT
javac -d "$classes_dir" "$repo_root/tests/java_reference/ToadsRelationFactorsReferenceHarness.java"
java -cp "$classes_dir" ToadsRelationFactorsReferenceHarness > "$repo_root/tests/data/synthetic/toads_relation_factors/java_reference.json"
