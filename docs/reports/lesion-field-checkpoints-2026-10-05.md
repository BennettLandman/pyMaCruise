# LesionTOADS adaptive-field checkpoint integration — 2026-10-05

Implementation commit: `cd6d3f2cd76b0583bad67466a22c34085c8ffe24`.

The registered-loop checkpoint can now accept an optional channel-last
correction-field NIfTI on the registered T1/FLAIR RAS grid. It hashes the input,
applies the source-selected adaptive field refit during loop passes, and stores
the resulting fields and current polynomial degree in the checkpoint summary
and output manifest. A later loop segment restores both from the parent
checkpoint. Field correction remains opt-in, so existing workflow defaults
retain their previous behavior.

The restartable `plus-lesion-loop-run` workflow includes field options in its
input/settings signature and passes them to the registered-loop stage. Both
that command and `plus-lesion-registered-loop` expose field path, degree ceiling,
and initial degree controls. Output field images are `.nii.gz`, RAS, on the
registered modalities grid, and protected by the checkpoint's output hashes.

Focused command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py
```

Result: **11 passed, 0 failed, 0 skipped in 125.71s**. Tests exercise a
synthetic initial field, adaptive degree-one execution, persisted RAS fields,
degree/field continuation into a second segment, workflow resume/repair, and
CLI argument routing.

Full-suite evidence for this same implementation commit is in
[the full regression report](full-pytest-2026-10-05-lesion-field-checkpoints.md):
**1,201 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,205.29s**;
peak sampled RSS was **10.986 GB** under the 24 GB watchdog. This is synthetic
Python integration and selected-method component evidence, not composite Java
equivalence or human-MRI validation. Selected base-TOADS Java line accounting
is unchanged by this LesionTOADS-only slice.
