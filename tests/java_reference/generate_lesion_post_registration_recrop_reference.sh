#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
class_dir=${TMPDIR:-/tmp}/pymacruise-lesion-recrops
mkdir -p "$class_dir"
javac -d "$class_dir" \
  "$repo_root/tests/java_reference/LesionPostRegistrationRecropReferenceHarness.java"
java -cp "$class_dir" LesionPostRegistrationRecropReferenceHarness \
  > "$repo_root/tests/data/synthetic/lesion_preprocessing_crop_rigid/post_registration_recrop_java_reference.json"
