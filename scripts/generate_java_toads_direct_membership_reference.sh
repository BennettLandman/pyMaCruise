#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-direct-membership-classes"
output="$repo_root/tests/data/synthetic/toads_direct_membership/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsDirectMembershipReferenceHarness.java"
java -cp "$classes" ToadsDirectMembershipReferenceHarness > "$output"
