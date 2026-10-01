#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
fixture="$repo/tests/data/synthetic/lesion_preprocessing_crop_rigid"
tmp="$repo/build/java-reference/lesion-preprocessing-crop-rigid"
mkdir -p "$tmp" "$fixture"
javac -d "$tmp" "$repo/tests/java_reference/LesionPreprocessingCropRigidReferenceHarness.java"
java -cp "$tmp" LesionPreprocessingCropRigidReferenceHarness > "$fixture/java_reference.json"
shasum -a 256 "$repo/tests/java_reference/LesionPreprocessingCropRigidReferenceHarness.java" "$fixture/java_reference.json"
