#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-channel-weights-classes"
output="$repo_root/tests/data/synthetic/lesion_channel_weights/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" \
  "$repo_root/tests/java_reference/LesionChannelWeightsReferenceHarness.java"
java -cp "$classes" LesionChannelWeightsReferenceHarness > "$output"
