#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
ref="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-p18czp-actual-next-loop-classes"
mkdir -p "$classes"
javac -d "$classes" \
  "$ref/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java" \
  "$ref/ToadsP18czpActualNextLoopMembershipReferenceHarness.java"
java -cp "$classes" ToadsP18czpActualNextLoopMembershipReferenceHarness "$@"
