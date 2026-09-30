# Selected LesionTOADS brainstem distance fixture

This synthetic slice compares the exact atlas-name selection in
`LesionToadSegmentation.computeBstemDistance()` followed by its byte-valued
`fastDistance(boolean[][][], byte, int, int, int)` helper. It does not invoke
the full LesionTOADS workflow or validate lesion segmentation.

`seeded_labels.nii.gz` is an asymmetric `(7, 8, 9)` label volume with one
border `Brainstem`, interior `InterVentricular-WM` and `FornixWM` voxels, and
three alias/control names. It is compared with caps 1 and 4. The separate
`aliases_only_labels.nii.gz` volume has no exact-name seeds and exercises the
empty-mask sentinel case. Both files have RAS axis codes, anisotropic voxel
sizes, and nonzero origins. The Java source arrays use `(x, y, z)` indexing;
the flattened JSON arrays preserve that order.

Only these exact names are selected: `Brainstem`, `InterVentricular-WM`, and
`FornixWM`. The lowercase, missing-hyphen, and extra-hyphen variants are
negative controls, matching Java `String.equals`. Distance uses six-connected
voxel steps, ignores physical spacing, is capped, and uses `-1` for unreachable
voxels. The selected kernel's `cap=1` extra propagation behavior is retained.

Reference provenance:

- Source file: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`
- Source SHA-256: `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`
- Excerpt methods: `computeBstemDistance()` and byte `fastDistance(...)`
- Standalone reference harness: `tests/java_reference/LesionToadsBrainstemDistanceReferenceHarness.java`
- Harness SHA-256: `0de0717502f1601244a0d498122bf3b0d15a0df2368d1c76952835300411d97f`
- Generated `java_reference.json` SHA-256: `b93533e9d4029fd6f15e477b8cf75375927582003a1741bfe1ac237b21d89ab4`

Regenerate the Java reference and NIfTI phantoms from the repository root:

```sh
mkdir -p /tmp/pymacruise-brainstem-reference
javac -d /tmp/pymacruise-brainstem-reference \
  tests/java_reference/LesionToadsBrainstemDistanceReferenceHarness.java
java -cp /tmp/pymacruise-brainstem-reference \
  LesionToadsBrainstemDistanceReferenceHarness \
  > tests/data/synthetic/lesion_brainstem_distance/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_brainstem_distance/generate_fixtures.py
```

The NIfTI generator writes only the two named files in this directory and
does not delete or modify data elsewhere. The reference commands create a
compiler-output directory under `/tmp` and overwrite only the named reference
JSON and fixture files inside this fixture directory; they contain no cleanup
or deletion commands.
