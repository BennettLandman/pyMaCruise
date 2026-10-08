# TOADS pyramid nonfinite-trial policy — 2026-10-08

## Failure found in full-suite evaluation

The full suite on tree `1d562a0` completed with **1,499 passed, 2 failed, 0 skipped, and 11 expected RAS warnings** in **1,786.29s**. Peak sampled RSS was **11.597 GB** under the 16 GB watchdog; it did not terminate the run. The failures were:

- `test_registered_startup_runs_all_three_source_scales_at_four_levels` still asserted the helper-transcription result (NaN coarse energy). With the actual VandyCRUISE method semantics, this fixture has a finite final rejected-trial energy and a finite transform handoff.
- `test_toads_loop_cli_runs_supported_multilevel_pyramid_after_nan_trial` exposed a real source edge case. The Java `Numerics.floor(NaN)` and comparison-based `bounded` path can turn NaN coordinates into a finite edge interpolation. On this fixture, the selected Java pyramid accepts a NaN transform at scale 4, then throws `ArrayIndexOutOfBoundsException` at scale 2 while sampling derivatives. The initial Python correction evaluated these trials the same way but stopped earlier with a finite-transform validation error.

## Implemented behavior

The Python pyramid now rejects nonfinite trial transforms by default, retaining the last finite transform for finer levels. The actual Java interpolation behavior remains available with `allow_nonfinite_trials=True` in the Python API and `--source-faithful-nonfinite-registration-trials` in `toads-loop`; this opt-in can reproduce the historical failure on degenerate input. This keeps ordinary workflow runs usable while making the source behavior explicit for controlled comparison.

## Focused evidence

- `pytest -q tests/test_toads_four_level_registration.py tests/test_toads_loop_cli.py tests/test_toads_registration_sampling.py tests/test_toads_alignment_startup.py` — **34 passed in 65.62s** on implementation commit `d42ee16`.
- The direct actual-method Java comparison still sets `allow_nonfinite_trials=True` and matches final transform, lambda, energy, coarse rejected-trial trace, and scale-2 trial counts.
- A full-suite rerun after this policy change is pending. Until it passes, P18ck remains Started and this is not a full-suite acceptance claim.

The source comparison remains component-level synthetic evidence, not whole-workflow Java/MATLAB equivalence or real-MRI validation.
