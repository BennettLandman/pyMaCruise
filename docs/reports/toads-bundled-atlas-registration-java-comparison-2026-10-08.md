# Bundled-atlas TOADS startup registration comparison

Date: 2026-10-08  
Scope: selected VandyCRUISE actual rigid pyramid method using the bundled 11-class atlas and synthetic brain-resolution T1/mask inputs.

## Reproduction

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_toads_actual_atlas_registration.py \
    --max-iterations 50 \
    --output-dir build/toads-actual-atlas-registration-50iter-2026-10-08
```

The script loads the repository atlas, prepares the committed synthetic `(172,216,164)` T1/mask startup, and compares its registration against `ToadsPyramidActualRegistrationHarness`, which calls the selected `MultipleToadDeformableAtlas` methods directly. The Java and Python runs both use a 50-iteration maximum at each of the three source scales.

## Result

- Exit code **0**; selected Java and Python final six-parameter float32 transforms match bit-for-bit.
- Eleven atlas classes; startup membership grid `(167,207,160,11)` and shape-prior grid `(172,216,164,11)`.
- At scales 8, 4, and 2, Python naturally stops after 12, 7, and 9 optimizer iterations, respectively; no 50-iteration cap is reached. All three stop with `energy_tolerance` and finite transforms.
- Python's accepted updates and total line-search trials match Java at each scale: 11/46, 6/30, and 9/35. The terminal transform is bit-exact.
- Peak sampled RSS: **1.140 GB** under a 16 GB watchdog.
- Reproducibility outputs: `comparison.json`, `java-harness.log`, and raw harness arrays under the requested output directory. The report records T1, brain-mask, and atlas-descriptor SHA-256 hashes.

This is synthetic actual-atlas registration evidence with natural convergence. It does not validate full TOADS evolution, the two-pass MaCRUISE composite, real MRI, or atlas-to-subject biological alignment. It adds no selected Java source lines to the coarse coverage estimate.
