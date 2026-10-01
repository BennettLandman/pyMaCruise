#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
BUILD="$ROOT/build/java-reference/toads-inhomogeneity-basis"
OUT="$ROOT/tests/data/synthetic/toads_inhomogeneity_basis/java_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" "$ROOT/tests/java_reference/ToadsInhomogeneityBasisReferenceHarness.java"
java -cp "$BUILD" ToadsInhomogeneityBasisReferenceHarness "$ROOT/tests/data/synthetic/toads_inhomogeneity_basis/grid.nii.gz" > "$OUT"
