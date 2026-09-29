#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/lesion-toads-membership-integration-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_membership_integration/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$reference_root/LesionToadsMembershipIntegrationReferenceHarness.java"
java -cp "$classes" LesionToadsMembershipIntegrationReferenceHarness > "$output"
