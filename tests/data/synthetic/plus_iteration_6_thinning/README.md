# P41z synthetic iteration-6 thinning inputs

This fixture freezes the completed iteration-5 registered state used as input
to the next LesionTOADS `initThinningLabels` + `propagateThinning` pass. It
includes the matching segmentation, skeleton, memberships, availability,
ordering, and the 13 registered shape-prior maps. All volumes are `.nii.gz`
on the same RAS grid; `provenance.json` records each file hash and the source
checkpoint hashes.

The shape maps were rebuilt from the iteration-5 rigid transform and support
bounds using the class-ordered shape maps in the committed atlas descriptor,
the registered modality affine for image spacing, and the descriptor atlas
affine for template spacing. The resulting 71×89×67×13 float32 array is
bit-exact to the saved `registered_recrop/transformed_shapes_13.nii.gz` output.
The fixture records both voxel-array hashes and the original saved-shape-file
hash. The 17 membership channels are ordered as 13 class maps, lesion, black
hole, outlier, and relation; the selected thinning methods consume the first
15 channels.

The focused pytest compares the Python result with the committed Java
`Lesion13ClassNextThinningReferenceHarness`. It requires a JDK and uses only
the committed fixture, Java excerpt, heap helper, and topology LUT. The byte
comparison covers initialization volumes and boundary queue, every pop and
insert event, candidate/topology decisions, remaining queue, segmentation,
skeleton, availability, ordering, and final class volumes. This is selected
method equivalence for one synthetic checkpoint, not validation of the whole
LesionTOADS iteration or workflow.

```sh
PYTHONPATH=src .venv/bin/pytest tests/test_plus_iteration6_thinning_java_reference.py -q
```

To reconstruct the input fixture from the ignored workflow checkpoint, which
is not needed to run the pytest:

```sh
PYTHONPATH=src .venv/bin/python scripts/materialize_plus_iteration6_thinning_fixture.py
```

After a successful comparison, the Python post-thinning state can be captured
separately for the next statistics comparison by setting
`PYMACRUISE_CAPTURE_ITERATION6_THINNING_DIR` to an output directory. The
capture writes RAS segmentation, skeleton, availability, and ordering NIfTIs
plus class volumes and trace counts in `result.json`; the test never changes
the committed input fixture.
