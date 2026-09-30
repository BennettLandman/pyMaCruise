#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
vandy_root="$repo_root/../toads-cruise-scm-2026-09-26/VandyCRUISE/src"
seg_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
atlas_source="$vandy_root/edu/jhmi/rad/medic/methods/MultipleToadAtlas.java"
power_source="$vandy_root/edu/jhmi/rad/medic/utilities/PowerTable.java"
descriptor="$repo_root/data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt"
[ "$(shasum -a 256 "$seg_source" | awk '{print $1}')" = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884" ] || { echo "Unexpected MultipleToadSegmentation source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$atlas_source" | awk '{print $1}')" = "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db" ] || { echo "Unexpected MultipleToadAtlas source hash" >&2; exit 1; }
[ "$(shasum -a 256 "$descriptor" | awk '{print $1}')" = "096eb94f9c814f47535e308b8b0d63eaacbbfc646efc1dd06c8ecd8f0892f006" ] || { echo "Unexpected TOADS2012 descriptor hash" >&2; exit 1; }
data_dir="$repo_root/tests/data/synthetic/toads_atlas_type_membership"
classes_dir="$repo_root/build/java-reference/toads-atlas-type-membership"
mkdir -p "$classes_dir"
javac -d "$classes_dir" "$power_source" "$repo_root/tests/java_reference/ToadsAtlasTypeMembershipReferenceHarness.java"
java -Xmx768m -cp "$classes_dir" ToadsAtlasTypeMembershipReferenceHarness \
    "$data_dir/images.nii.gz" "$data_dir/intensity_field.nii.gz" "$data_dir/psf_bias.nii.gz" \
    "$data_dir/initial_memberships.nii.gz" "$data_dir/transformed_shapes.nii.gz" \
    "$data_dir/closest_labels.nii.gz" "$data_dir/folded_ribbon.nii.gz" "$descriptor" \
    > "$data_dir/java_reference.json"
