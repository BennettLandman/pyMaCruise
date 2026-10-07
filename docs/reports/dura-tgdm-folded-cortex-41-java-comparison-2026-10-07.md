# Dura TGDM on a 41³ folded cortex — selected Java comparison

## Scope

This extends the selected VandyCRUISE `DynamicPressureTGDM` comparison to a
41×41×41 analytic folded-cortex phantom. The selected Java harness runs the
VandyCRUISE TGDM/profile sources with `DistanceFieldyh`; it does not execute
the complete `MedicAlgorithmRemoveDura09` wrapper or its topology-correction
stages. It is synthetic stage evidence, not full Dura/workflow equivalence or
real-MRI validation. Inputs and Java outputs are committed as RAS `.nii.gz`
files in [`tests/data/synthetic/dura_tgdm_folded_cortex_41`](../../tests/data/synthetic/dura_tgdm_folded_cortex_41),
with a SHA-256 manifest and source/parameter provenance.

The five-pass configuration is the selected Dura TGDM profile:
`membership_level=0.4`, `alpha=0.5`, `alignment_force=1.0`,
`curvature_force=0.5`, `pressure_force=0.5`, 150 updates. The Java harness
reports 30 updates per pass and zero landmine restarts.

## Results

All five Python/Java post-reinitialization pass fields have identical signs
(Dice 1.0). Their maximum absolute-value differences are 0.00788, 0.00930,
0.00537, 0.01453, and 0.01049; mean absolute differences are `1.79e-6`,
`7.81e-6`, `9.29e-6`, `2.32e-5`, and `1.38e-5`. The final level-set array has
zero sign disagreements, Dice **1.0**, maximum absolute difference
**0.01049**, and mean absolute difference **1.38e-5**.

The final meshes have 1,572 vertices each. Symmetric mean nearest-surface
distance is **0.0000816 mm** (Python→Java `0.0000838 mm`, Java→Python
`0.0000795 mm`); Hausdorff distance is **0.00350 mm**. Both are within the
project gates (Dice ≥0.97, symmetric mean distance ≤0.01 mm, Hausdorff ≤1 mm).

The largest intermediate difference is the fourth pre-reinitialization field:
maximum absolute error **0.34674**, but no sign changes. The error is not
hidden in the final surface metrics: the report and saved pass fields preserve
its location in the algorithm sequence. It is a known source-numerical
variation requiring stage-specific review; final labels and geometry still
meet the stated gates on this fixture.

## Reproduction

From the repository root, with a JDK:

```sh
JAVA_TOOL_OPTIONS=-Xmx2g \
DURA_TGDM_NX=41 DURA_TGDM_NY=41 DURA_TGDM_NZ=41 \
DURA_TGDM_FOLDED=true \
DURA_TGDM_OUTPUT=build/java-reference/dura-folded-41/java_reference.json \
scripts/generate_java_dura_tgdm_full_evolution_reference.sh

./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_dura_folded_41_java.py \
  --import-java-json build/java-reference/dura-folded-41/java_reference.json

./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/pytest -q tests/test_dura_folded_cortex_41.py
```

The focused regression passed **1 test in 5.25 s**, peak sampled RSS **0.149
GB**. The full-crop comparison and fixture conversion passed under the same 8 GB
watchdog. No full repository suite was run for this slice.

## Remaining Dura gap

This validates only TGDM evolution on one larger folded synthetic case. The
full Dura wrapper still needs a selected-source composite comparison that
covers the two topology-correction calls, WM+GM union, initial distance fields,
TGDM, and final mask/output transfer in source order, then integration with the
brain-resolution base workflow. The stage result must not be described as full
Dura equivalence.
