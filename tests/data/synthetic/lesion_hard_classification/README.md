# LesionTOADS hard-classification method fixture

This fixture validates only the selected `LesionToadSegmentation.generateClassification()` voxelwise map operation. It maps the highest positive ordinary membership through `templateLabel`; the lesion label replaces it only when lesion membership is strictly greater. It does not validate any LesionTOADS iteration or workflow output wrapper.

`memberships_and_lesion.nii.gz` is a synthetic 4-D RAS NIfTI of shape `(10, 1, 1, 4)`, with unit voxel spacing and a zero origin. The final axis stores three ordinary membership channels followed by the lesion membership. The ten voxels exercise a strict lesion win, an ordinary tie, all-zero scores, negative-only scores, NaN handling, positive-infinity ties, equality between lesion and ordinary maxima, lesion-only positive support, ordinary class ties, and a NaN lesion score. Values are synthetic and contain no human imaging data.

The Java harness emits the expected signed-byte label sequence in `java_reference.json`. It includes a negative signed-byte atlas label (`-106`) to exercise label preservation as well as labels that differ from class indices. Its source method is pinned to `LesionToadSegmentation.java` SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; the copied selected-method excerpt SHA-256 is `885d221e7414fc76fd1000c9dca36a1725643b4ba5aa18c686567611c0d0c10a`, and the reference JSON SHA-256 is `4a9deec7593852552c1a1f6c5baa9b0120c4e07f26008b2aaad0fc75c5c37064`. The compressed synthetic NIfTI SHA-256 is `bd351d39da63b3431eff206d3b8c4dc9fe9ae3ab5aa267a78d435907a0e0c0bc`.

Regenerate the reference from the checked-in harness and the input image with:

```sh
mkdir -p /tmp/pymacruise-lesion-hard-classification
javac -d /tmp/pymacruise-lesion-hard-classification \
  tests/java_reference/LesionToadsHardClassificationReferenceHarness.java
java -cp /tmp/pymacruise-lesion-hard-classification \
  LesionToadsHardClassificationReferenceHarness \
  > tests/data/synthetic/lesion_hard_classification/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_hard_classification/generate_fixtures.py
```

The generator writes only `memberships_and_lesion.nii.gz` in this fixture directory.
