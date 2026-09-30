# LesionTOADS objective reduction fixture

`objective_inputs_ras.nii.gz` is a deterministic float32 4-D synthetic image
with shape `5 x 6 x 7 x 11` and a positive RAS affine. Its channels are three
already-normalized class memberships, three class energy terms, the source
denominator, lesion membership, black-hole membership, lesion energy, and
black-hole energy. The generator includes positive and zero class energies,
zero/negative/NaN denominator cases, a positive-infinity term, and an
infinity-times-zero NaN case.

`LesionToadsObjectiveReductionReferenceHarness.java` transcribes the selected
active source reduction order for interior voxels: visit x/y/z, add class
energy products in ascending class order only when `den > 0`, then evaluate
the grouped lesion-plus-black-hole energy expression and add it unconditionally.
The reference stores the Java float32 running total after each voxel as well
as the final total. IEEE `NaN`/infinity values are intentionally retained.

The harness is driven by a text export of the exact NIfTI float values so Java
does not need a NIfTI dependency. Rebuild both fixture and reference with:

```sh
.venv/bin/python tests/data/synthetic/lesion_objective_reduction/generate_fixture.py
.venv/bin/python tests/data/synthetic/lesion_objective_reduction/generate_reference.py
```

The pinned selected VandyCRUISE source is
`LesionToadSegmentation.java`, SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.

This helper only reduces supplied objective terms with supplied normalized
maps and denominator values. It does not calculate data, smoothing, or prior
energy terms; estimate statistics; normalize memberships; or run the full
LesionTOADS workflow.
