#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-registration-volume-classes"
output="$repo_root/tests/data/synthetic/lesion_registration_volume/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionRegistrationVolumeAccumulationReferenceHarness.java"
java -cp "$classes" LesionRegistrationVolumeAccumulationReferenceHarness > "$output"
