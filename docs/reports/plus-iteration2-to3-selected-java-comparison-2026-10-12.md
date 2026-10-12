# PLUS synthetic iteration 2→3 selected Java comparison — 2026-10-12

## Result

The continuation from the one-pass checkpoint (iteration 2) through iteration 3 matches the compared selected Java methods and the saved Python iteration-3 checkpoint. This is a held-fixed synthetic stage comparison, not full Java workflow equivalence or human-data validation.

The compared behavior matched exactly: registration shape sampling (60,880,560 float32 values), 13-class thinning serialization (136,921,848 bytes), post-thinning statistics, next GM membership fields and NaN counts, registration transform/rotation/bounds, post-registration growth (1,099,864 events and all checkpoint arrays), and relation refresh (127 relations, 14 diffusion iterations). The checkpointed step changes 119,918 segmentation voxels relative to iteration 2; this state transition is reproduced by the selected Java method comparisons. The membership update contains **25,552,556 NaNs in both Java and Python**, including the source-faithful NaN objective; the energy ratio is NaN in both. Finite float32 fields are compared bit-for-bit, and NaN values are compared by classification plus exact count rather than payload/sign bits.

Despite this selected-method match, the subsequent surface audit worsens: GM support drops from 82.98% at iteration 2 to 73.59% at iteration 3, with failed connected-genus-zero, nesting, and production-thickness checks. Thus, the observed degradation is consistent with the source-selected state evolution on this synthetic input, not evidence by itself of a mismatch in the compared stages. The broader input-to-surface composition and synthetic fixture’s anatomical suitability remain open.

The first comparator attempt exposed report-only defects: strict JSON rejected source-faithful NaN scalars, and the Python NaN tally omitted the objective while Java includes it. The comparator now writes nonfinite scalars as explicit strings and counts the objective NaN; no algorithm values were changed. The complete comparison then exited successfully under a 16 GB RSS watchdog with **7.924 GB peak sampled RSS**.

## Reproduction

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/plus-labeled-folded-deterministic-numba-onepass-2026-10-12/lesion_workflow/registered-loop/lesion-registered-loop.json \
  --output-dir build/plus-labeled-folded-java-iteration2-to3-2026-10-12 \
  --expected-next-checkpoint \
  build/plus-labeled-folded-deterministic-numba-twopass-2026-10-12/lesion_workflow/registered-loop/lesion-registered-loop.json
```

The committed machine-readable comparison is [plus-iteration2-to3-selected-java-comparison-2026-10-12.json](plus-iteration2-to3-selected-java-comparison-2026-10-12.json). Both checkpoints are local ignored `build/` artifacts from the deterministic 65³ folded synthetic case. Their manifests hash the inputs and outputs; the compared starting manifest SHA-256 is recorded in the JSON. The run used source-faithful zero-variance NaN behavior and RAS outputs.
