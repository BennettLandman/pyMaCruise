# Selected LesionTOADS ventricle distance fixture

This synthetic slice compares the exact name selection in the historically
misspelled `LesionToadSegmentation.computeVentirclsDistance()` method and its
byte-valued `fastDistance(boolean[][][], byte, int, int, int)` helper. It does
not invoke the full LesionTOADS workflow or validate lesion segmentation.

`seeded_labels.nii.gz` is an asymmetric `(7, 8, 9)` integer label volume. It
contains all four exact selected atlas names, with two border seeds and two
interior seeds, plus case, punctuation, and unrelated negative controls. The
same input checks caps 1 and 4. `nonmatching_labels.nii.gz` contains only
negative controls and checks the empty-mask sentinel result with cap 3. Both
files have true RAS axis codes, anisotropic spacing, and a nonzero origin. The
Java source array and JSON reference use `(x, y, z)` indexing and C-order
flattening. The header asserts RAS; the numerical comparison uses source voxel
indices, as this selected distance kernel ignores physical spacing.

Only these exact spellings select voxels: `Ventricles`, `LeftVentricle`,
`RightVentricle`, and `ThirdVentricle`. Names are not case- or punctuation-
normalized, matching Java `String.equals`. The distance output is a signed
byte array using six-connected voxel steps, `-1` for unreachable voxels, and
the Java cap behavior (including one propagation layer beyond a cap of 1).

Reference provenance:

- Source file: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`
- Source SHA-256: `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`
- Excerpt methods: `computeVentirclsDistance()` and byte `fastDistance(...)`
- Standalone reference harness:
  `tests/java_reference/LesionToadsVentricleDistanceReferenceHarness.java`
- Harness SHA-256: `aa4f8a7ccddd3644bae6b9a675efc08da30aa20269d36f927d3517b34d5da45c`
- Generated `java_reference.json` SHA-256:
  `d886c02988386465684ac014b1c444d3a609c129a0cc950972004f794104dd62`

Regenerate the Java reference and NIfTI inputs from the repository root:

```sh
mkdir -p /tmp/pymacruise-ventricle-reference
javac -d /tmp/pymacruise-ventricle-reference \
  tests/java_reference/LesionToadsVentricleDistanceReferenceHarness.java
java -cp /tmp/pymacruise-ventricle-reference \
  LesionToadsVentricleDistanceReferenceHarness \
  > tests/data/synthetic/lesion_ventricle_distance/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_ventricle_distance/generate_fixtures.py
```

These commands create compiler files under `/tmp` and overwrite only the
reference JSON and two fixture volumes inside this directory.
