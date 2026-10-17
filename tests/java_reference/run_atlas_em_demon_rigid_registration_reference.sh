#!/bin/sh
set -eu

if [ "$#" -ne 7 ]; then
  echo "usage: $0 nx ny nz iterations levels single|pyramid arrays-directory" >&2
  exit 2
fi

reference_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "$reference_root/../.." && pwd)
source_root=${ATLAS_EM_SOURCE_ROOT:-"$(CDPATH= cd -- "$project_root/../mipav-src" && pwd)/src"}
toads_root=${TOADS_SOURCE_ROOT:-"$(CDPATH= cd -- "$project_root/../toads-cruise-scm-2026-09-26/toads-cruise" && pwd)/src"}
classes="$project_root/build/java-reference/atlas-em-demon-rigid"

check_hash() {
  actual=$(shasum -a 256 "$1" | awk '{print $1}')
  if [ "$actual" != "$2" ]; then
    echo "unexpected selected source hash for $1: $actual" >&2
    exit 1
  fi
}

check_hash "$toads_root/edu/jhmi/rad/medic/methods/DemonToadDeformableAtlas.java" 96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a
check_hash "$toads_root/edu/jhmi/rad/medic/methods/ParametricTransform.java" 0071ac49c6f7b612345534120aa7fbe561da316169f96c60fbe6d06ce0d37216
check_hash "$toads_root/edu/jhmi/rad/medic/structures/RotationMatrix.java" 144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0
check_hash "$toads_root/edu/jhmi/rad/medic/libraries/ImageFunctions.java" 0b017e75be715c6e5842fa1788a331cf9cc64416b170c7d2844cf75818f6e98a
check_hash "$toads_root/edu/jhmi/rad/medic/utilities/Numerics.java" f3baf2cc372c164cb0194c0396416cbb627828b50df226407f4817d6fc851fe2
check_hash "$toads_root/edu/jhmi/rad/medic/methods/DemonsToadsAtlasWarping.java" 7037749ffe502b9476f56cd90a5d6ad7ff3da9140a79835d802660e8e38efe28
mkdir -p "$classes"
javac -J-Xmx1g -sourcepath "$toads_root:$reference_root/stubs:$reference_root/upstream" -d "$classes" \
  "$reference_root/stubs/gov/nih/mipav/view/MipavUtil.java" \
  "$reference_root/stubs/gov/nih/mipav/model/structures/jama/HeadlessPlaceholder.java" \
  "$reference_root/stubs/edu/jhu/ece/iacl/utility/ArrayUtil.java" \
  "$reference_root/upstream/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$reference_root/AtlasEMDemonRigidRegistrationReferenceHarness.java" \
  "$toads_root/edu/jhmi/rad/medic/methods/DemonToadDeformableAtlas.java"
java -Xmx1g -cp "$classes" AtlasEMDemonRigidRegistrationReferenceHarness "$@"
