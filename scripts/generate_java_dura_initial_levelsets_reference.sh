#!/bin/sh
set -eu

if [ "${JIST_TEST_JAR:-}" = "" ]; then
  echo "Set JIST_TEST_JAR to the supplied original JIST-Test jar." >&2
  exit 2
fi

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/dura-initial-levelsets-classes"
output="$repo_root/tests/data/synthetic/dura_initial_levelsets/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -cp "$JIST_TEST_JAR" -d "$classes" \
  "$repo_root/tests/java_reference/DuraInitialLevelSetsJarHarness.java"
java -cp "$classes:$JIST_TEST_JAR" DuraInitialLevelSetsJarHarness "$output"
