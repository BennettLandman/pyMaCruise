# Terminal LesionTOADS registered-loop checkpoint

## Implemented behavior

P41m now persists a source convergence or iteration-limit branch instead of
rejecting every stopped iteration. Terminal checkpoints write the final
template-label classification as signed int8 and the lesion-subtype map as
signed int16, both as RAS `.nii.gz`, and record stop flags/reason, iteration,
classification counts, and output hashes. Continuing checkpoints retain
their skeleton, availability, ordering, and relation outputs and remain
restartable. A terminal checkpoint can be verified/resumed as a completed
result, but it cannot be hydrated as an input to another segmentation
iteration.

P41o accepts a `--loop-convergence-threshold` independently from P41k's
`--max-distance`, reports when the overall run produced a terminal result, and
returns the terminal classification/subtype paths. If a bounded segment does
not stop, it still returns the continuing checkpoint state.

## Validation

- A synthetic P41m test exercises both terminal outcomes: an iteration limit at
  iteration 2 and convergence at iteration 2 using an explicit synthetic
  threshold. Both persisted arrays match `LesionToadsFinalClassification`
  voxelwise and retain their required data types and RAS affines.
- The terminal checkpoint test verifies manifest hashes, stop metadata,
  rejection as a continuation input, and cache reuse on resume. Existing
  continuing checkpoint tests still verify resume and corrupt-output repair.
- P41o synthetic end-to-end now runs P41i → P41k → P41m through an
  iteration-limit terminal output; top-level manifest and returned terminal
  paths are checked. CLI argument routing and help were updated for the
  separate convergence threshold.
- Focused `tests/test_plus_lesion_registered_branch.py`: **10 passed** in
  103.65s; the workflow CLI routing selector passed separately (**1 passed**
  in 0.68s).
- Full suite on implementation commit
  `0ace0ffc330bf733de0a4ebebfb0d9cccb2b61eb`: **1,192 passed, 0 failed, 0
  skipped, 11 expected RAS warnings in 1,014.58s**. Peak sampled RSS was
  **14.168 GB** under the 24 GB watchdog. See the [full-suite report](full-pytest-2026-10-05-plus-terminal-lesion.md).

This is synthetic Python workflow evidence. The final-classification adapter
has source-level tests, but these checkpoint tests are not a complete Java
plus-run comparison or real-MRI validation.
