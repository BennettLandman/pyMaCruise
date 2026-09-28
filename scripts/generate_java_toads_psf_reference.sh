#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/toads-psf-bias-classes"
output="$repo_root/tests/data/synthetic/toads_psf_bias/java_reference.json"

mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/ToadsPsfBiasReferenceHarness.java"
java -cp "$classes" ToadsPsfBiasReferenceHarness > "$output"
