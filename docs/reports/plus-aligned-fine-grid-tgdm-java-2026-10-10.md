# Fine-grid synthetic plus TGDM Java comparison — 2026-10-10

## Scope

This is a held-fixed comparison of the selected `MaGenericTGDM` wrapper with
the saved Python TGDM outputs from the iteration-6 checkpoint of the labeled
81×101×77 (~2 mm) synthetic plus case. The comparison is on the workflow's
73×91×63 RAS crop. Hash-checked crop inputs were used to
reconstruct the tissue, ACE, and GVF fields; the saved initial signed-distance
level set was reused exactly. TOADS/LesionTOADS and all upstream stages remain
held fixed. The Java process peaked at **0.538 GB sampled RSS** under an 8 GB
watchdog. It is synthetic source-comparison evidence, not full-workflow
equivalence or real-MRI validation.

## Results

The saved Python and Java initial level sets match exactly. The three final
surface masks each differ at only three voxels, with Dice ≥0.99993. However,
the agreed Hausdorff limit is **1 mm**, and the measured distances exceed it:

| Surface | Sign-different voxels | Dice at zero | Mean surface distance | Hausdorff | Gate |
|---|---:|---:|---:|---:|---|
| Central | 3 | 0.999935 | 0.000371 mm | 1.345 mm | Fail |
| Inner | 3 | 0.999934 | 0.000375 mm | 1.374 mm | Fail |
| Outer | 3 | 0.999935 | 0.000510 mm | 1.347 mm | Fail |

The same three coordinates already differ in sign in the first inner surface,
so the final-surface Hausdorff miss is inherited from that early TGDM stage.
This is a small but source-visible topology-sensitive discrepancy; it cannot be
reported as passing the Java gate or hidden by the Dice result. The independent
ribbon-quality audit still fails its tissue-support and thickness criteria.

## Reproduction

The held-fixed inputs and Java products are generated under ignored `build/`
directories. Recreate the inputs and compile/run the selected Java harness:

```sh
mkdir -p build/plus-vandy-atlas-aligned-fine/java-reference/classes

./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --context-manifest build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/synthetic-plus-checkpoint-surface-composition.json \
  --output-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --use-saved-tgdm-outputs

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
  precomputed_initial_level_set \
  > build/plus-vandy-atlas-aligned-fine/java-inputs/java-watchdog.log 2>&1

./.venv/bin/python scripts/compare_plus_tgdm_java.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-tgdm-java-2026-10-10.json
```

The complete measurements and hashes are in the adjacent machine-readable
[JSON report](plus-aligned-fine-grid-tgdm-java-2026-10-10.json).
