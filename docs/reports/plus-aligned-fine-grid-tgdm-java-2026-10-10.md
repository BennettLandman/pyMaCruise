# Fine-grid synthetic plus TGDM Java comparison — 2026-10-10

## Scope

This is a held-fixed comparison of the selected `MaGenericTGDM` wrapper with
the saved Python TGDM outputs from the iteration-6 checkpoint of the labeled
81×101×77 (~2 mm) synthetic plus case. The comparison is on the workflow's
73×91×63 RAS crop. Hash-checked crop inputs were used to reconstruct the
tissue, ACE, and GVF fields. Java and Python start from the same held-fixed
initial-membership array and independently construct the signed-distance
seed. The resulting seed fields have identical signs and a maximum absolute
difference of `9.54e-7`. TOADS/LesionTOADS and all other upstream stages remain
held fixed. The Java process peaked at **0.604 GB sampled RSS** under an 8 GB
watchdog. This is synthetic component source-comparison evidence, not
full-workflow equivalence or real-MRI validation.

## Results

The Java and Python initial level sets have no sign differences. The three
final surface masks each differ at only three voxels, with Dice ≥0.99993. Production distances below use the source-matched JIST connectivity-18 extractor used by PyMaCRUISE. Marching-cubes values are retained in the JSON as diagnostics only. The agreed Hausdorff limit is **1 mm**, and the initial-inner and later surfaces exceed it:

| Surface | Sign-different voxels | Dice at zero | Mean surface distance | Hausdorff | Gate |
|---|---:|---:|---:|---:|---|
| Smoothed initial | 1 | 0.999986 | 0.000184 mm | 0.313 mm | Pass |
| Initial inner | 3 | 0.999934 | 0.000334 mm | 1.411 mm | Fail |
| Central | 3 | 0.999935 | 0.000343 mm | 1.414 mm | Fail |
| Inner | 3 | 0.999934 | 0.000335 mm | 1.412 mm | Fail |
| Outer | 3 | 0.999935 | 0.000406 mm | 1.415 mm | Fail |

The smoothed initial surface differs in sign at one voxel
`(29, 37, 26)` (Java `0.000082`, Python `-0.000670`), while its measured
surface distances still pass. The first inner surface is the first stage to
miss the source gate; its three sign disagreements persist in all later
surfaces. This is a small but source-visible topology-sensitive discrepancy;
it cannot be reported as passing the Java gate or hidden by the Dice result.
The independent ribbon-quality audit still fails its tissue-support and
thickness criteria.

As a diagnostic, Python was rerun from the exact Java smoothed level set through
the initial inner, central, recomputed inner, and outer stages. All four surfaces pass the Java thresholds with zero sign differences and Dice 1.0 using production JIST surfaces:

| Surface | Mean surface distance | Hausdorff |
|---|---:|---:|
| Initial inner | 0.00000062 mm | 0.001588 mm |
| Central | 0.00000479 mm | 0.003462 mm |
| Inner | 0.00000058 mm | 0.000453 mm |
| Outer | 0.00000437 mm | 0.002890 mm |

The diagnostic uses source-matched JIST mesh extraction for its primary gate; marching-cubes values remain secondary diagnostics. It deliberately substitutes a Java intermediate and is not an
end-to-end Python-vs-Java result. It localizes the amplified mismatch to the
Python smoothing output passed into the inner stage. The explicit Java
multiply/divide order for curvature magnitude is now encoded in the vectorized
kernel; a full replay still reproduced all six saved Python products exactly,
and the held-fixed comparison metrics did not change. The focused TGDM suite
passed **26 tests** after this arithmetic alignment. The downstream diagnostic
is reproducible with the committed
[comparison script](../../scripts/compare_plus_surfaces_from_java_smoothing.py)
and its [machine-readable report](plus-aligned-fine-grid-java-smoothing-downstream-diagnostic-2026-10-10.json).

The comparison harness now saves `smoothed_initial` and the comparator measures
it and `initial_inner`, making this failure location reproducible. An
exploratory run with float32 rounding enabled for signed-distance updates did
not change the sign-disagreement counts, so that precision option is not
recorded as a correction.

## Reproduction

The held-fixed inputs and Java products are generated under ignored `build/`
directories. Recreate the inputs and compile/run the selected Java harness:

```sh
mkdir -p build/plus-vandy-atlas-aligned-fine/java-reference/classes

./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --context-manifest build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/synthetic-plus-checkpoint-surface-composition.json \
  --output-dir build/plus-vandy-atlas-aligned-fine/java-inputs

javac -d build/plus-vandy-atlas-aligned-fine/java-reference/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java

./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g \
  -cp build/plus-vandy-atlas-aligned-fine/java-reference/classes \
  PlusTGDMReferenceHarness \
  build/plus-vandy-atlas-aligned-fine/java-inputs/arrays.bin \
  build/plus-vandy-atlas-aligned-fine/java-inputs/java \
  > build/plus-vandy-atlas-aligned-fine/java-inputs/java-watchdog.log 2>&1

./.venv/bin/python scripts/compare_plus_tgdm_java.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-tgdm-java-2026-10-10.json

./.venv/bin/python scripts/compare_plus_surfaces_from_java_smoothing.py \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-java-smoothing-downstream-diagnostic-2026-10-10.json
```

The complete measurements and hashes are in the adjacent machine-readable
[JSON report](plus-aligned-fine-grid-tgdm-java-2026-10-10.json).
