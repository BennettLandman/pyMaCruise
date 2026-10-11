# P41z fine-grid Dura reference case

This committed fixture preserves the exact registered LesionTOADS state after
iteration 6 from the synthetic P41z run. It contains 15 registered membership
channels (atlas classes 0–12, lesion membership 13, black-hole membership 14),
the matching segmentation, and the four Python Dura outputs. All images are
`.nii.gz` on the checkpoint's RAS grid. `provenance.json` records hashes of the
original ignored checkpoint files and hashes of each committed image.

The selected Java method is
`LesionToadSegmentation.generateDuraRemovalOutputs`; its source SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
The Java harness reads these committed arrays directly. The comparison test
uses the committed harness and does not need the original Java source checkout.
The fixture is not a
real MRI, and exact Dura output agreement does not establish equivalence of
the upstream LesionTOADS iterations or the complete MaCRUISEplus workflow.

Rebuild the fixture only when the original ignored P41z checkpoint is present:

```sh
PYTHONPATH=src .venv/bin/python scripts/materialize_fine_grid_lesion_dura_fixture.py
```

Run the comparison from the committed fixture:

```sh
.venv/bin/python -m pytest tests/test_plus_fine_grid_dura_java_reference.py -q
```
