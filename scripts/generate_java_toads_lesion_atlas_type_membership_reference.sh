#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
vandy_root="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src"
seg_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
atlas_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadAtlas.java"
power_source="$vandy_root/edu/jhmi/rad/medic/utilities/PowerTable.java"
descriptor="$repo_root/data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
[ "$(shasum -a 256 "$seg_source" | awk '{print $1}')" = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884" ] || { echo "Unexpected MultipleToadSegmentation source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$atlas_source" | awk '{print $1}')" = "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db" ] || { echo "Unexpected MultipleToadAtlas source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$power_source" | awk '{print $1}')" = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678" ] || { echo "Unexpected PowerTable source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$descriptor" | awk '{print $1}')" = "e42901c860901e0c562e01e8800eda6ac43c528ea4c7bdf27218ad8daab870d7" ] || { echo "Unexpected lesion TOADS descriptor hash" >&2; exit 1; }
data_dir="$repo_root/tests/data/synthetic/toads_lesion_atlas_type_membership"
classes_dir="$repo_root/build/java-reference/toads-lesion-atlas-type-membership"
mkdir -p "$classes_dir"
javac -d "$classes_dir" "$power_source" "$repo_root/tests/java_reference/ToadsLesionAtlasTypeMembershipReferenceHarness.java"
java -Xmx768m -cp "$classes_dir" ToadsLesionAtlasTypeMembershipReferenceHarness \
    "$data_dir/images.nii.gz" "$data_dir/intensity_field.nii.gz" "$data_dir/psf_bias.nii.gz" \
    "$data_dir/initial_memberships.nii.gz" "$data_dir/transformed_shapes.nii.gz" \
    "$data_dir/closest_labels.nii.gz" "$data_dir/segmentation.nii.gz" \
    "$data_dir/initial_lesion_membership.nii.gz" "$data_dir/folded_ribbon.nii.gz" "$descriptor" \
    > "$data_dir/java_reference.json"
