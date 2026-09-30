#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-default-registration-handoff-classes"
output="$repo_root/tests/data/synthetic/toads_default_registration_handoff/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsDefaultRegistrationHandoffReferenceHarness.java"
java -cp "$classes" ToadsDefaultRegistrationHandoffReferenceHarness > "$output"
