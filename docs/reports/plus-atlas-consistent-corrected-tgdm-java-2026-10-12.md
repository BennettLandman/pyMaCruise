# Corrected synthetic PLUS TGDM held-fixed Java comparison — 2026-10-12

## Scope

This is a **synthetic, held-fixed TGDM comparison**, not a whole-workflow Java comparison or independent anatomical validation. It compares the selected `MaGenericTGDM` plus-wrapper sequence against the saved Python PLUS TGDM stage on the corrected 172×216×164 atlas-aligned case. The source signed-distance seed and saved TOADS/LesionTOADS products are reused; boundary, tissue, ACE, and GVF fields are reconstructed from the checked workflow crop inputs. This isolates TGDM evolution from upstream changes.

The harness uses the already-saved TGDM products when assembling held-fixed inputs. Its `smoothed_initial` result is replayed separately because the public stage manifest does not persist that intermediate. Java ran with the explicit `precomputed_initial_level_set` mode; omitting that mode would incorrectly reinterpret the signed level set as a membership map.

## Results

All five agreed surface tolerances pass: zero-level Dice ≥0.97, symmetric mean surface distance ≤0.01 mm, and sampled symmetric Hausdorff distance ≤1 mm. `initial_level_set` is bit-exact. `smoothed_initial`, `initial_inner`, `central`, and `inner` have Dice 1.0 and sub-0.00002 mm maximum sampled surface distances. `outer` is the limiting output: Dice **0.999943**, mean surface distance **0.001837 mm**, and Hausdorff **0.937940 mm**. It has 121 sign-changed voxels, but remains inside the agreed surface-position limits. Peak Java RSS was **1.197 GB** under an 8 GB watchdog.

| Output | Bit-exact | Changed voxels | Mean absolute voxel error | Maximum absolute voxel error | Dice at zero | Mean surface distance (mm) | Hausdorff (mm) | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| initial level set | yes | 0 | 0 | 0 | — | — | — | exact |
| smoothed initial | no | 283,850 | 3.38×10⁻⁸ | 1.79×10⁻⁵ | 1.000000 | 3.54×10⁻⁸ | 1.50×10⁻⁵ | pass |
| initial inner | no | 108,809 | 1.04×10⁻⁸ | 2.83×10⁻⁶ | 1.000000 | 8.55×10⁻⁹ | 1.41×10⁻⁵ | pass |
| central | no | 140,232 | 1.14×10⁻⁸ | 4.41×10⁻⁶ | 1.000000 | 8.39×10⁻⁹ | 1.26×10⁻⁵ | pass |
| inner | no | 79,276 | 5.12×10⁻⁹ | 1.87×10⁻⁵ | 1.000000 | 4.28×10⁻⁹ | 1.34×10⁻⁵ | pass |
| outer | no | 2,314,773 | 0.001460 | 0.974887 | 0.999943 | 0.001837 | 0.937940 | pass |

These findings narrow the cause of the terminal ribbon failures: the selected Python TGDM stage is within the agreed Java-reference tolerances on the same inputs. The remaining differences must be investigated in upstream TOADS/LesionTOADS, Dura transfer, the input tissue construction, or the anatomy/acceptance assumptions. The outer surface has less than 0.063 mm of Hausdorff margin, so preserve this comparison as a regression case.

## Reproduction

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 \
  --report-interval-seconds 60 -- .venv/bin/python \
  scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/plus-atlas-consistent-corrected-full-2026-10-12/plus-surface \
  --output-dir build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12 \
  --use-saved-tgdm-outputs --allow-python-stage-differences

javac -d build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java

PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 2 \
  --report-interval-seconds 30 -- /usr/bin/java -Xms256m -Xmx6g \
  -cp build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12/classes \
  PlusTGDMReferenceHarness \
  build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12/arrays.bin \
  build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12/java \
  precomputed_initial_level_set

PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/compare_plus_tgdm_java.py \
  --workflow-root build/plus-atlas-consistent-corrected-full-2026-10-12/plus-surface \
  --reference-dir build/plus-atlas-consistent-corrected-java-tgdm-saved-2026-10-12 \
  --report docs/reports/plus-atlas-consistent-corrected-tgdm-java-2026-10-12.json
```

The adjacent JSON contains per-voxel errors, sign changes, surface metrics, source-mode metadata, and run signatures. This pass validates a TGDM component only; it does not resolve upstream label support, TOADS convergence, or whole-workflow parity.
