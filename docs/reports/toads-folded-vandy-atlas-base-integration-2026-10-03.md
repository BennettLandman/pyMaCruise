# Folded Vandy atlas synthetic base integration — 2026-10-03

## Purpose and evidence boundary

This experiment advances the synthetic base-workflow checkpoint with a
folded-anatomy fixture derived from the VandyCRUISE atlas already bundled in
the repository. The fixture downsamples the source atlas maps into a compact
grid, synthesizes T1 intensities from the descriptor's T1_MPRAGE class priors,
and supplies a matching brain mask and MaSeg label map. The stored fixture is
compressed NIfTI in RAS orientation; it contains no human imaging data.

This is Python synthetic integration evidence. It is not Java whole-loop
equivalence, clinical validation, or evidence that full-resolution atlas
spacing behaves identically. The selected TOADS stages in the restartable base
run are deliberately capped at two iterations, so their `converged` values
must remain false.

## Reproduce the fixture

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py
./.venv/bin/pytest -q tests/test_toads_downsampled_vandy_atlas_fixture.py
```

The generator resamples the atlas image arrays in their shared voxel-index
grid. This is necessary because the historical topology NIfTI header says
1 mm while its source descriptor and aligned shape maps use 0.9375 mm. A
header-based interpolation would move the topology away from the shape maps;
the generator preserves the source's raw shared array layout and writes the
descriptor's resampled spacing to all generated assets.

## Bounded convergence/runtime probe

The source-atlas loop command used the selected wrapper values
(`relation_limit=3`, descriptor optimized factors `[0.8, 1.0, 1.0]`,
`prior_scale=0.2`, `prior_factor=0.5`, and growing enabled):

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --fixture tests/data/synthetic/toads_downsampled_vandy_atlas \
  --atlas tests/data/synthetic/toads_downsampled_vandy_atlas/atlas.txt \
  --output build/downsampled-vandy-atlas-10iter.nii.gz \
  --maximum-iterations 10 --registration-iterations 1 \
  --loop-registration-iterations 1 --padding 2 --relation-limit 3
```

Result: alignment took 4.582 seconds and the 10-iteration loop took
187.337 seconds. The first and last objectives were 100.346458 and 58.3192406.
The energy ratio fell to 0.00312617 at iteration 9, then rose to 0.0180569 at
iteration 10. The output was a 36×45×34 uint8 RAS NIfTI with all 11 source
labels and SHA-256
`cdfc966e6106e2e66a5817cd0f3f8d097f251d5c2205c10263c8dfbb1c0f026f`.
It did not reach the `1e-4` threshold within ten iterations. This is a
markedly different trajectory from the tiny procedural atlas probe, which
grew to extreme energies; no conclusion about full 50-iteration convergence
is drawn until that run completes.

The Python exact membership update took about 9.5–10.0 seconds per iteration
on this 36×45×34 phantom. Registration and branch work added about 8 seconds
per continuing cycle. Runtime remains a practical blocker to repeated full
source-cap experiments; the stopping threshold has not been loosened.

## Restartable synthetic base run

The first attempted run was interrupted after creating an empty TOADS stage
directory. That exposed a restart bug: recovery picked a fresh directory but
still checked outputs under the old planned attempt path. The runner now
preserves orphan attempts, selects the next unused attempt number, and remaps
expected outputs into the selected directory. The source atlas materializer
also accepts the selected Vandy structure names as well as the synthetic
`Csf1/Gm1/Wm1` aliases.

The resumed command completed the full six-stage base route:

```bash
./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_downsampled_vandy_atlas/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_downsampled_vandy_atlas/brain-mask.nii.gz \
  --maseg tests/data/synthetic/toads_downsampled_vandy_atlas/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_downsampled_vandy_atlas/atlas.txt \
  --output-dir build/downsampled-vandy-base-run \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1 --resume
```

Manifest `build/downsampled-vandy-base-run/run.json` is `complete`, records
all six stages as complete, and correctly records
`stage-run-complete-but-toads-did-not-meet-energy-threshold`. Both TOADS
stages ran two iterations and reached their iteration cap, with terminal
energy ratios 0.244500503 and 0.247157827. The workflow wrote the three
required nested level sets, cropped tissue and ACE/GVF intermediates, refined
segmentation, provenance, stage logs, and a self-contained derived atlas.
The second invocation with `--resume` verified all six stages and skipped
their computations. NIfTI outputs are RAS and shape-checked against the
native T1 grid, except the explicitly recorded CRUISE crop-grid products.

## Validation and next integration milestone

- `./.venv/bin/pytest -q tests/test_base_workflow_restartable.py tests/test_toads_downsampled_vandy_atlas_fixture.py`
  — **3 passed in 74.06s** on implementation commit `7fb4886`. This includes
  source-atlas descriptor handling and restart from a pre-existing orphan
  attempt directory.
- The 10-iteration loop trace is a diagnostic command run on a synthetic
  atlas-derived phantom, not a pytest reference and not Java equivalence.
- The whole base run passed all stage checks and a live resume verification;
  its two TOADS stages are capped and therefore not a converged workflow.
- No full pytest suite or real-data run was performed for this change.

Next: finish the source-cap synthetic convergence/runtime experiment, determine
whether vectorizing or compiling the sequential membership sweep can reduce
runtime while preserving Java float/order behavior, and then use the same
fixture for a converged restartable base run before advancing plus.
