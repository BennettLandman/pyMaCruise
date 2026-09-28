#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-powertable-classes"
output="$repo_root/tests/data/synthetic/toads_powertable/java_reference.json"
source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$source" "$repo_root/tests/java_reference/ToadsPowerTableReferenceHarness.java"
java -cp "$classes" ToadsPowerTableReferenceHarness > "$output"
