#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-optimized-membership-pass-classes"
output_dir="$repo_root/tests/data/synthetic/toads_optimized_membership_pass"
power_source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"

mkdir -p "$classes" "$output_dir"
javac -d "$classes" "$power_source" "$repo_root/tests/java_reference/ToadsOptimizedMembershipPassReferenceHarness.java"
java -cp "$classes" ToadsOptimizedMembershipPassReferenceHarness false 0.0 > "$output_dir/java_reference_field_off_psf_off.json"
java -cp "$classes" ToadsOptimizedMembershipPassReferenceHarness false 0.2 > "$output_dir/java_reference_field_off_psf_on.json"
java -cp "$classes" ToadsOptimizedMembershipPassReferenceHarness true 0.0 > "$output_dir/java_reference_field_on_psf_off.json"
java -cp "$classes" ToadsOptimizedMembershipPassReferenceHarness true 0.2 > "$output_dir/java_reference_field_on_psf_on.json"
