#!/bin/sh
set -eu

if [ "$#" -ne 23 ]; then
  echo "usage: $0 mem.bin shapes.bin nix niy niz ntx nty ntz classes rix riy riz rtx rty rtz levels iterations tx ty tz rx ry rz" >&2
  exit 2
fi

reference_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "$reference_root/../.." && pwd)
source_root=${VANDYCRUISE_SOURCE_ROOT:-"$(CDPATH= cd -- "$project_root/.." && pwd)/toads-cruise-scm-2026-09-26/VandyCRUISE/src"}
classes="$project_root/build/java-reference/toads-pyramid-actual-registration"
check_hash() {
  actual=$(shasum -a 256 "$1" | awk '{print $1}')
  if [ "$actual" != "$2" ]; then
    echo "unexpected selected Java source hash for $1: $actual" >&2
    exit 1
  fi
}
check_hash "$source_root/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java" e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2
check_hash "$source_root/edu/jhmi/rad/medic/methods/ParametricTransform.java" 22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac
check_hash "$source_root/edu/jhmi/rad/medic/structures/RotationMatrix.java" 144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0
check_hash "$source_root/edu/jhmi/rad/medic/libraries/ImageFunctions.java" be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892
mkdir -p "$classes"
javac -sourcepath "$source_root:$reference_root/stubs" -d "$classes" \
  "$reference_root/stubs/gov/nih/mipav/view/MipavUtil.java" \
  "$reference_root/stubs/gov/nih/mipav/model/structures/jama/HeadlessPlaceholder.java" \
  "$reference_root/upstream/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java" \
  "$reference_root/ToadsPyramidActualRegistrationHarness.java" \
  "$source_root/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java"
java -cp "$classes" ToadsPyramidActualRegistrationHarness "$@"
