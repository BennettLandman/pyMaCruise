#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
classes="$repo_root/build/java-reference/lesion-toads-membership-conversion-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_membership_conversion/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$upstream" "$reference_root/LesionToadsMembershipConversionReferenceHarness.java"
java -cp "$classes" LesionToadsMembershipConversionReferenceHarness > "$output"
