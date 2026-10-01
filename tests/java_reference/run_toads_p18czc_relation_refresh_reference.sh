#!/bin/sh
set -eu

if [ "$#" -ne 8 ]; then
  echo "usage: $0 segmentation-i8.bin memberships-f32le.bin boundary-i8.bin nx ny nz classes limit" >&2
  exit 2
fi

reference_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "$reference_root/../.." && pwd)
classes="$project_root/build/java-reference/toads-p18czc-relation-refresh"
mkdir -p "$classes"
javac -d "$classes" \
  "$reference_root/ToadsRelationListReferenceHarness.java" \
  "$reference_root/ToadsRelationDiffusionReferenceHarness.java" \
  "$reference_root/ToadsBoundaryDistancesReferenceHarness.java" \
  "$reference_root/ToadsClosestLabelsReferenceHarness.java" \
  "$reference_root/ToadsP18czcRelationRefreshReferenceHarness.java"
java -cp "$classes" ToadsP18czcRelationRefreshReferenceHarness "$@"
