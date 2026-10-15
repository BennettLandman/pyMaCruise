# P41z synthetic registered-loop continuation: iterations 90–99 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS continuation reached the configured
99-iteration limit. The final energy ratio was `0.13026626408100128`, far
above the source-selected `1e-4` natural-convergence threshold. The actual
limit-99 checkpoint records `stop_reason=iteration_limit`,
`iteration_limit_stop=true`, and `convergence_stop=false`. This is a clear
nonconvergence result for this selected synthetic state, not a claim that the
algorithm is unportable or that a complete original PLUS run has the same
trajectory.

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 90 | 0.1310552061 | No |
| 91 | 0.1299555749 | No |
| 92 | 0.1299000680 | No |
| 93 | 0.1303516030 | No |
| 94 | 0.1304168999 | No |
| 95 | 0.1310453415 | No |
| 96 | 0.1311262399 | No |
| 97 | 0.1307514310 | No |
| 98 | 0.1307056844 | No |
| 99 | 0.1302662641 | Iteration cap |

The 90–94 segment peaked at **1.200 GB** sampled RSS under a 4 GB limit; the
95–98 segment peaked at **1.141 GB**; iteration 99 peaked at **0.794 GB**.
All saved checkpoints verify their hashes and contain finite RAS NIfTIs: 9/9
and 8/8 for continuing checkpoints, 7/7 and 6/6 for the terminal product.

## Selected-method comparison at the limit

The terminal output contains no skeleton, availability, ordering, or relation
map, so the existing transition comparator cannot use it as a full
registered-loop checkpoint. To compare the 98→99 update, the same iteration
was rerun from the identical iteration-98 checkpoint with an iteration-100
ceiling, exposing the internal continuation state at iteration 99. This is a
method-level probe, not the actual selected cap-99 run. Both cap-99 and probe
runs produce the exact same iteration-99 energy ratio
`0.13026626408100128`.

The selected Java comparison of the 98→99 probe is exact for 5,503,849 shape
samples, 267,472 thinning events, statistics and membership/energy fields,
registration and bounds, 254,382 growth events, and 130 relation rows. Peak
comparator RSS was **0.958 GB** under 4 GB. This is selected-method parity on
matched checkpoint data; it does not compare the complete Java application
or independently validate Java's iteration-limit wrapper branch.

P41z remains an unfinished full-workflow validation row. A terminal surface
composition and explicit ribbon-quality audit were subsequently run from the
hash-verified iteration-99 state; those measurements are recorded below as
diagnostics, separate from source-equivalence evidence.

## Iteration-99 source-profile surface and ribbon audit

The actual terminal cap-99 checkpoint was composed through the selected Dura,
native-grid transfer, MaSeg separation, and CRUISE surface pipeline. It used
the committed labeled 81×101×77 RAS fixture and the source-selected TGDM
profile recorded in the composition manifest: initial/inner/central/outer
iterations **6/4/7/7**, WM iso `0.5`, GM iso `0.4`, central pressure `1.5`,
and the selected curvature/external-force settings. This is a composition from
saved checkpoints, not one uninterrupted PLUS invocation and not an
independent Java surface-evolution comparison. The loop's selected 98→99
methods and cap-control branch have separate Java evidence above.

The terminal composition completed at 0.634 GB peak sampled RSS under an 8 GB
watchdog. The explicit ribbon audit completed at 0.260 GB under a 4 GB
watchdog. It reports RAS geometry, 0 nesting violations, and watertight
manifolds with no boundary or nonmanifold edges. All four requested synthetic
MaSeg label values are present. The measured anatomy/geometry checks are:

- GM-between-surfaces support **93.03%**, below the diagnostic 97% reference;
  WM-inside-inner **99.42%** and CSF-outside-outer **99.27%**.
- Meshes have **2 connected components** each and Euler characteristic
  **−16**, so the connected genus-zero diagnostic fails.
- Level-set subtraction thickness mean is **4.841 mm** across 71,420 crossing
  voxels. Closest-surface thickness median/p05/p95 is
  **6.424/0.232/12.320 mm**; the median and tails miss the audit's nominal
  5 mm profile ranges.
- The audit JSON reports `passed=false` for tissue support, topology, and
  production thickness, with nesting and level-set mean checks passing.

These are diagnostic anatomy results, not criteria for declaring the port
incorrect; source agreement on matched inputs remains the acceptance test.
The audit records loop ratio `0.13026626408100128` and
`stop_reason=iteration_limit`, consistent with the separate nonconvergence
evidence. Machine output is under ignored
`build/p41z-case-normalization/iteration99-terminal-surface-audit/`;
composition signature `5dddf9bed372e20e238bd52390c7b2358289678217575bd20742740d08be825f`.

Reproduce from repository root:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/p41z-case-normalization/iteration99-terminal/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/p41z-case-normalization/iteration99-terminal-surface-audit

PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 1 --report-interval-seconds 15 -- \
  .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/p41z-case-normalization/iteration99-terminal-surface-audit/synthetic-plus-terminal-surface-composition.json \
  --report build/p41z-case-normalization/iteration99-terminal-surface-audit/ribbon-quality.json \
  --maseg-profile synthetic
```

## Same-input selected Java TGDM comparison

The TGDM harness was run on the exact cropped membership, tissue, ACE, and
GVF arrays reconstructed from that surface composition; it consumed the
hash-verified Python initial level set as its signed-distance seed. This is a
held-fixed TGDM comparison. Dura, native transfer, MaSeg separation, and the
loop state remain Python-generated inputs, so it is not full PLUS source
equivalence.

| Surface | Sign changes | Dice at zero | Symmetric mean distance | Sampled Hausdorff | Gates |
|---|---:|---:|---:|---:|---|
| Initial-inner | 0 | 1.000000 | 7.18×10⁻⁹ mm | 0.0000077 mm | Pass |
| Central | 0 | 1.000000 | 7.21×10⁻⁹ mm | 0.0000076 mm | Pass |
| Delivered inner | 0 | 1.000000 | 8.81×10⁻⁹ mm | 0.0000111 mm | Pass |
| Outer | 7 | 0.999970 | 0.001368 mm | 0.570943 mm | Pass |

The outer field is not voxel-identical (mean absolute field error
`0.000628`, maximum `0.321583`), but its zero-level surface remains within the
agreed Dice `≥0.97`, symmetric mean `≤0.01 mm`, and Hausdorff `≤1 mm` gates.
The 73×91×63 RAS comparison peaked at **0.667 GB** Java RSS under an 8 GB
watchdog. Comparison JSON and serialized arrays remain in ignored
`build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/`.

Reproduce from repository root (the third Java harness argument
`precomputed_initial_level_set` is required for this saved seed):

```sh
PYTHONPATH=src .venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/p41z-case-normalization/iteration99-terminal-surface-audit/plus-surfaces \
  --output-dir build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison \
  --context-manifest build/p41z-case-normalization/iteration99-terminal-surface-audit/synthetic-plus-terminal-surface-composition.json \
  --use-saved-tgdm-outputs

mkdir -p build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/java \
  build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/classes
javac -d build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 1 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g \
  -cp build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/classes \
  PlusTGDMReferenceHarness \
  build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/arrays.bin \
  build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/java \
  precomputed_initial_level_set \
  > build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/java-watchdog.log 2>&1
PYTHONPATH=src .venv/bin/python scripts/compare_plus_tgdm_java.py \
  --workflow-root build/p41z-case-normalization/iteration99-terminal-surface-audit/plus-surfaces \
  --reference-dir build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison \
  --report build/p41z-case-normalization/iteration99-terminal-surface-audit/tgdm-java-comparison/comparison.json
```

## Direct comparison to the saved iteration-99 terminal product

The comparator now accepts the iteration-99 cap manifest as the expected
next checkpoint and explicitly lists products that the terminal manifest
does not save. Replaying from the hash-verified iteration-98 checkpoint
matched Java exactly for shape sampling (5,503,849 values), next thinning
(267,472 queue events), statistics, membership and energy, registration
transform/bounds, growth (254,382 events), and relations (130 rows, 14
diffusion passes). The saved terminal segmentation matches the Java growth
result bit-for-bit. The terminal manifest has no skeleton, availability,
ordering, class-volume image, relation map, or relation-membership array, so
these checkpoint-output comparisons are reported as unavailable. This does
not compare the complete Java application or loop wrapper. A separate Java
control-branch harness evaluated the source-selected wrapper decision using
the actual iteration-98 and iteration-99 energies at threshold `1e-4`, with
iteration and cap both set to 99 and narrow-band mode disabled. The Java
float32 ratio bits (`3e056485`) equal Python's; Java reports
`convergence_stop=false`, `iteration_limit_stop=true`, `stop=true`,
`run_growing=true`, and `refresh_relations=false`, matching the Python cap
manifest and the selected 98→99 growth comparison. Harness source SHA256:
`dceb95289dea92bb97665fdb96b7bd71ad28dc2f75a3b4dfa7f4968cf444deb8`;
embedded source excerpt SHA256: `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.
This validates that control branch only, not orchestration or application-level
behavior. Focused comparator regressions: **2 passed**; combined with the Dura
comparator regression, **3 passed**.

The reproduction command compares the iteration-98 checkpoint directly to
the terminal cap-99 manifest:

```sh
.venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iterations95-98/lesion-registered-loop.json \
  --expected-next-checkpoint build/p41z-case-normalization/iteration99-terminal/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration99-terminal-java-comparison
```

## Reproduction

Continue the iteration-89 checkpoint to iteration 94 with
`pymacruise plus-lesion-registered-loop`, the bundled 12-object atlas
descriptor, `--iterations 99`, `--thinning-limit 1000`, and
`--max-passes 5`, under a 4 GB RSS watchdog. Resume from iteration 94 and use
`--max-passes 4` to iteration 98, then one final pass to 99. For the
selected-method probe only, repeat the final pass from the iteration-98
manifest with `--iterations 100`; compare its internal 99 state using
`scripts/compare_p41z_registered_thinning_java.py`. Artifacts remain under
ignored `build/p41z-case-normalization/`.
