#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/lesion-toads-fielded-data-energy-membership-classes"
output="$repo_root/tests/data/synthetic/lesion_toads_fielded_data_energy_membership/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$reference_root/LesionToadsDataEnergyMembershipReferenceHarness.java"
java -cp "$classes" LesionToadsDataEnergyMembershipReferenceHarness field > "$output"
