# Selected LesionTOADS gray-matter distance fixture

This synthetic slice compares the exact atlas-name selection in
`LesionToadSegmentation.computeGMDistance()` and its byte-valued
`fastDistance` overload. It does not exercise the complete LesionTOADS
segmentation workflow.

`seeded_atlas_indices.nii.gz` is a true-RAS `(9, 10, 11)` atlas-index volume.
It includes one voxel for each of the 13 accepted names, five case/punctuation
and unrelated controls, and multiple seeds including boundary seeds. The
same volume is checked with distance caps 1 and 4. `near_miss_atlas_indices.nii.gz`
contains only negative-control names and exercises the empty-mask result.
Both files have positive RAS axes, anisotropic spacing, and a nonzero origin;
distance remains six-connected voxel steps and ignores physical spacing.

Accepted names are exactly `Cerebrum-GM`, `CerebralGM`, `CorticalGM`,
`Caudate`, `GM`, `Thalamus`, `Putamen`, `Sulcal-CSF`, `SulcalCSF`,
`Subcortical-GM`, `SubcorticalGM`, `BasalGanglia`, and `GlobusPallidus`.
Near misses are not normalized. The byte map uses `-1` for unreached voxels
and preserves the source byte-kernel cap-1 extra-layer behavior.

Reference provenance:

- Source: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`
- Source SHA-256: `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`
- Excerpt methods: `computeGMDistance()` and byte `fastDistance(...)`
- Harness: `tests/java_reference/LesionToadsGmDistanceReferenceHarness.java`
- Harness SHA-256: `be18a8292350505187e34f669d4963b47f106359a0f513dd8202626e7e3768e3`
- Reference JSON SHA-256: `08e3aeabaec1588821f3e7b63e0c50117cdc39924f581804c1258bc57e6c559d`

Regenerate the Java reference and phantom from the repository root:

```sh
mkdir -p /tmp/pymacruise-gm-reference
javac -d /tmp/pymacruise-gm-reference tests/java_reference/LesionToadsGmDistanceReferenceHarness.java
java -cp /tmp/pymacruise-gm-reference LesionToadsGmDistanceReferenceHarness \
  > tests/data/synthetic/lesion_gm_distance/java_reference.json
PYTHONPATH=src .venv/bin/python tests/data/synthetic/lesion_gm_distance/generate_fixtures.py
```

The commands write only the named Java reference, NIfTI files, and compiler
output under `/tmp`; they do not remove files.
