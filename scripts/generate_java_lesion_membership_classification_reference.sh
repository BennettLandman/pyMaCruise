#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$repo/build/java-reference/lesion-membership-classification"
data="$repo/tests/data/synthetic/lesion_membership_classification"
mkdir -p "$build"
mkdir -p "$data"
javac -d "$build" "$repo/tests/java_reference/LesionMembershipClassificationReferenceHarness.java"
java -cp "$build" LesionMembershipClassificationReferenceHarness \
  > "$data/java_reference.json"
