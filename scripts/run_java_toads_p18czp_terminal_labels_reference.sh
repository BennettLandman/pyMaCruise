#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
ref="$repo_root/tests/java_reference/ToadsP18czpTerminalLabelsReferenceHarness.java"
classes="$repo_root/build/java-reference/toads-p18czp-terminal-labels-classes"
mkdir -p "$classes"
javac -d "$classes" "$ref"
java -cp "$classes" ToadsP18czpTerminalLabelsReferenceHarness "$@"
