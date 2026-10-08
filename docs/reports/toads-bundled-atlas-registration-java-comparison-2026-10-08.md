# Bundled-atlas TOADS startup registration comparison

Date: 2026-10-08  
Scope: selected VandyCRUISE actual rigid pyramid method using the bundled 11-class atlas and synthetic brain-resolution T1/mask inputs.

## Reproduction

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_toads_actual_atlas_registration.py \
    --output-dir build/toads-actual-atlas-registration-2026-10-08
```

The script loads the repository atlas, prepares the committed synthetic `(172,216,164)` T1/mask startup, performs one registration iteration per source pyramid level, and compares its final transform with `ToadsPyramidActualRegistrationHarness`. The harness calls the selected `MultipleToadDeformableAtlas` methods directly.

## Result

- Exit code **0**; selected Java and Python final six-parameter float32 transforms match bit-for-bit.
- Eleven atlas classes; startup membership grid `(167,207,160,11)` and shape-prior grid `(172,216,164,11)`.
- All three Python block scales (8, 4, 2) return finite transforms and accept one update each; Java/Python final transform bits match.
- Peak sampled RSS: **1.140 GB** under a 16 GB watchdog.
- Reproducibility outputs: `comparison.json`, `java-harness.log`, and raw harness arrays under the requested output directory.

This is a capped synthetic atlas-registration comparison. It validates the source mapping with the bundled production atlas, not 50-iteration convergence, full TOADS evolution, the two-pass MaCRUISE composite, real MRI, or atlas-to-subject biological alignment. It adds no selected Java source lines to the coarse coverage estimate.
