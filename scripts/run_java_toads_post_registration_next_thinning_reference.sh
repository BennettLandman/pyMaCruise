#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_root="$repo_root/tests/java_reference"
classes="$repo_root/build/java-reference/toads-post-registration-next-thinning-classes"
upstream="$reference_root/upstream/edu/jhmi/rad/medic/structures"

if [ "$#" -ne 5 ]; then
  echo "usage: $0 input.bin nx ny nz classes" >&2
  exit 2
fi

mkdir -p "$classes"
javac -d "$classes" \
  "$upstream/CriticalPointLUT.java" \
  "$upstream/BinaryHeap4D.java" \
  "$reference_root/ToadsHomeomorphicChangeReferenceHarness.java" \
  "$reference_root/ToadsPostRegistrationNextThinningReferenceHarness.java"
cp "$repo_root/src/pymacruise/data/critical186LUT.raw.gz" "$classes/critical186LUT.raw.gz"
java -cp "$classes" ToadsPostRegistrationNextThinningReferenceHarness "$@"
