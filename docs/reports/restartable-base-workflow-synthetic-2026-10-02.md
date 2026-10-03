# Restartable synthetic base workflow

## Runnable command

From the repository root, with ANTsPy installed, run:

```bash
./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --output-dir build/synthetic-base-run \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

Use `--resume` with the same inputs and settings to verify and reuse completed
stages. If a recorded output is missing or its SHA-256 changes, that stage and
dependent stages rerun in a new attempt directory. Earlier attempt files are
preserved for diagnosis.

## Implemented route and outputs

Commit `e27f0b6` adds `pymacruise base-run` and the synthetic MaSeg input. The
runner composes existing CLIs in this order:

1. ANTsPy N4 correction
2. Registered TOADS pass and Dura input products
3. Dura removal
4. Second registered TOADS pass
5. MaSeg-guided brainstem/cerebellum cleanup
6. Crop, topology correction, ACE, GVF, and nested TGDM
7. MaSeg segmentation refinement

`run.json` records input/atlas hashes, parameters and signature, stage commands,
output hashes, dependency signatures, attempt history, stage checks, and the
two TOADS stop summaries. Each stage log is also preserved. The expected final
products include native-T1-grid RAS `inner.nii.gz`, `central.nii.gz`,
`outer.nii.gz`, `SEG_refine.nii.gz`, and the segmentation-refinement rule
summary. The manifest includes intermediate `.nii.gz` images and provenance.

Focused validation:

```text
./.venv/bin/python -m pytest tests/test_base_workflow_restartable.py -q
1 passed in 88.62s
```

The test ran all seven stages, validated the native-grid finite outputs and
manifest, resumed a completed run with all stages verified, corrupted a saved
core level set, and resumed again. Verified upstream stages were skipped;
CRUISE core and dependent refinement reran in attempt 2. Scientific output
hashes matched the clean run. The N4 → two TOADS → Dura → CRUISE composition
and its source stop summaries are synthetic integration evidence, not Java
end-to-end equivalence. A later full repository suite includes this runner's
tests: **1,119 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,273.87s** on commit `963794e`; see the [TOADS optimization full-suite
checkpoint](full-pytest-2026-10-02-toads-runtime-optimization.md).

## Convergence status and boundary

Both TOADS stages in the test intentionally use the two-iteration compact
setting. The outputs are mechanically complete, but both reports say
`energy_threshold_reached=false` and `iteration_limit_reached=true`. The run
manifest sets `converged=false`; it must not be interpreted as a converged or
numerically accepted MaCRUISE segmentation. The selected Java stopping formula
is preserved. The separate [TOADS trace](synthetic-toads-loop-characterization-2026-10-02.md)
shows persistent energy oscillation on the current 13³ fixture and no
justification to relax the threshold.

No human MRI was used. Full source-equivalence and real-data validation remain
open.
