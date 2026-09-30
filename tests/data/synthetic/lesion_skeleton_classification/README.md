# LesionTOADS skeleton classification fixture

This fixture compares only `LesionToadSegmentation.updateClassificationFromSkeleton()`. The selected Java `EMPTY=-1` sentinel clears a classification voxel to zero; other skeleton bytes index the source template-label array. The method mutates the classification image and returns the count of changed bytes. This is a bounded adapter check, not a workflow equivalence claim; source call sites do not invoke this helper in the selected wrapper.

`initial_classification.nii.gz` and `skeleton.nii.gz` are synthetic RAS `(10, 1, 1)` int8 label volumes with unit spacing and zero origin. Cases include EMPTY, every available class index, signed template byte labels (including -128), unchanged and changed prior labels, and zero output. `java_reference.json` records the selected-method excerpt output and source provenance.

Provenance: selected `LesionToadSegmentation.java` SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; harness SHA-256 `8f98f9c1306fc09d83c382b07899fb1e96f6e93b8dbe115c1b6cc2e4c2fa7689`; reference JSON SHA-256 `2a6568450b2fe269d2c139940ed7b5038eef139e77520cd1a83f7d1a4bb0988a`; initial-classification and skeleton NIfTI SHA-256 values are `0856d24c97b6eb148c23cc7b5d716b97b3bda1145ab7adecbc1c57a94573262a` and `2d021028280513d3dc6eb004be0b74972398b6d5157d45a1816beabc791ff90a`, respectively.

Regenerate the synthetic inputs and reference from the in-repository harness with:

```sh
mkdir -p /tmp/pymacruise-lesion-skeleton-classification
javac -d /tmp/pymacruise-lesion-skeleton-classification \
  tests/java_reference/LesionToadsSkeletonClassificationReferenceHarness.java
java -cp /tmp/pymacruise-lesion-skeleton-classification \
  LesionToadsSkeletonClassificationReferenceHarness \
  > tests/data/synthetic/lesion_skeleton_classification/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_skeleton_classification/generate_fixtures.py
```

The generator writes only the two named NIfTI files under this fixture directory.
