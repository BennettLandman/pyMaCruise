#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
vandy_root="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src"
toads_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
power_source="$vandy_root/edu/jhmi/rad/medic/utilities/PowerTable.java"
expected_toads="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
expected_power="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
[ "$(shasum -a 256 "$toads_source" | awk '{print $1}')" = "$expected_toads" ] || { echo "Unexpected MultipleToadSegmentation source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$power_source" | awk '{print $1}')" = "$expected_power" ] || { echo "Unexpected PowerTable source hash" >&2; exit 1; }
data_dir="$repo_root/tests/data/synthetic/toads_default_atlas_outlier"
classes_dir="$repo_root/build/java-reference/toads-default-outlier"
mkdir -p "$classes_dir"
javac -d "$classes_dir" "$power_source" "$repo_root/tests/java_reference/ToadsDefaultOutlierMembershipReferenceHarness.java"
java -Xmx768m -cp "$classes_dir" ToadsDefaultOutlierMembershipReferenceHarness \
    "$data_dir/images.nii.gz" "$data_dir/field.nii.gz" "$data_dir/psf_bias.nii.gz" \
    "$data_dir/initial_memberships.nii.gz" "$data_dir/transformed_shapes.nii.gz" \
    "$data_dir/closest_labels.nii.gz" "$data_dir/folded_ribbon.nii.gz" \
    "$data_dir/initial_outlier_membership.nii.gz" > "$data_dir/java_reference.json"
