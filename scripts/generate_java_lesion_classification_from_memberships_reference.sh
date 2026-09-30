#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$repo/build/java-reference/lesion-classification-from-memberships"
data="$repo/tests/data/synthetic/lesion_classification_membership_update"
mkdir -p "$build" "$data"
javac -d "$build" "$repo/tests/java_reference/LesionClassificationFromMembershipsReferenceHarness.java"
java -cp "$build" LesionClassificationFromMembershipsReferenceHarness > "$data/java_reference.json"
