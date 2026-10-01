#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
UPSTREAM="$ROOT/tests/java_reference/upstream/edu/jhmi/rad/medic"
BUILD="$ROOT/build/java-reference/toads-lesion-approx-registration"
OUT="$ROOT/tests/data/synthetic/toads_lesion_approximate_registration/java_reference.json"
mkdir -p "$BUILD" "$(dirname -- "$OUT")"
javac -d "$BUILD" \
  "$UPSTREAM/utilities/MedicUtilPublic.java" \
  "$UPSTREAM/utilities/Numerics.java" \
  "$UPSTREAM/structures/RotationMatrix.java" \
  "$UPSTREAM/methods/ParametricTransform.java" \
  "$UPSTREAM/libraries/ImageFunctions.java" \
  "$ROOT/tests/java_reference/ToadsLesionApproximateRegistrationReferenceHarness.java"
java -cp "$BUILD" ToadsLesionApproximateRegistrationReferenceHarness \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/images.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/transformed_shapes.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/atlas_template_shapes.nii.gz" \
  "$ROOT/data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt" > "$OUT"
python3 -c 'import hashlib,json,sys; out,*paths=sys.argv[1:]; r=json.load(open(out)); names=("harness","numerics","images","transformed_shapes","template_shapes","descriptor"); r["p18czn_evidence_sha256"]={n:hashlib.sha256(open(p,"rb").read()).hexdigest() for n,p in zip(names,paths)}; open(out,"w").write(json.dumps(r,separators=(",",":"))+"\n")' \
  "$OUT" \
  "$ROOT/tests/java_reference/ToadsLesionApproximateRegistrationReferenceHarness.java" \
  "$UPSTREAM/utilities/Numerics.java" \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/images.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/transformed_shapes.nii.gz" \
  "$ROOT/tests/data/synthetic/toads_lesion_atlas_type_membership/atlas_template_shapes.nii.gz" \
  "$ROOT/data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
