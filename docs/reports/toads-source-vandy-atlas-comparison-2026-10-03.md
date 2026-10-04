# TOADS source-atlas comparison — 2026-10-03

## Scope

This is a synthetic Java-source/Python comparison on the committed
`(36,45,34)` Vandy-atlas-derived T1 fixture. The selected Java algorithm reads
headerless atlas payloads in their original descriptor index order. Python's
first comparison used canonical RAS arrays and diverged at crop/startup. The
source-order diagnostic then reproduced the Java crop and startup state. That
behavior is now used by the production TOADS registration path. This validates
the captured source slice only; it does not establish complete TOADS or
MaCRUISE equivalence, convergence on brain-resolution data, or real-data
performance.

## Atlas conversion and ordering

The conversion audit covers 35 bundled atlas images. The source-payload
round-trip check restored all 32 raw payloads exactly using their source
dimensions, datatype, and MIPAV voxel order. The three Analyze images were
voxel-exact after canonicalizing through their source header affine. Their
source-axis metadata still disagrees with the audit (`L,A,S` in the headers vs
`L,P,S` in the audit), so physical origin/direction claims for those examples
remain unverified. Those SPECTRE images are not used in this Vandy comparison.

The implementation now restores audited NIfTI payloads to Java index order
for TOADS registration. Standalone NIfTI atlases without a conversion audit
retain their native array order, using their affine axis codes. The T1 and
final outputs remain on the RAS input grid; the internal atlas index order
preserves the source algorithm's descriptor-based mapping.

## Matched startup and first-cycle results

Java used the original VandyCRUISE 11-structure raw descriptor and the
instrumented `VandyCoreSyntheticRunExport` harness. The synthetic T1 values
were exported in Fortran order to preserve MIPAV's x-fastest voxel convention.
Python used the committed converted NIfTI atlas through the production loader,
one registration level, one registration update per branch, padding 10, and
relation limit 3. Python's trace requires at least two loop iterations; the
captured first iteration is compared with Java's one-iteration run.

- Java and Python crop bounds match exactly: `[3,35] [3,42] [0,30]`; both use
  the `(53,60,51)` padded grid.
- The initial membership objective matches: Java `171.02353`, Python
  `171.02352905273438`.
- The source final rigid transform and Python post-refresh transform agree
  within `1.7e-7` across the six values (the Java diagnostic prints six
  significant digits).
- Transformed shape priors: mean absolute difference `1.0611e-8`, maximum
  `2.8611e-6`; 4,568 of 1,783,980 values exceed `1e-6`.
- Initial memberships: mean absolute difference `1.0744e-9`, maximum
  `1.9670e-6`; 66 of 1,783,980 values exceed `1e-6`.
- First exact membership update: mean absolute difference `5.6475e-10`,
  maximum `3.7849e-6`; 67 of 1,783,980 values exceed `1e-6`.
- Transformed constructor topology matches all 162,180 voxels exactly.
- First thinning labels match all 162,180 voxels exactly.
- Iteration 1 objective and energy ratio match: Java `50.30009 / 1.0909222`,
  Python `50.3000908 / 1.09092224`.

The reported values are selected Java-reference comparisons on one synthetic
input. They are not a Dice/Hausdorff/mean-surface-distance acceptance result
for the final MaCRUISE workflow. Python's final trace output is uint8 on the
input's `(36,45,34)` RAS grid. No real MRI was used.

## Evidence and remaining work

- Java source log: ignored diagnostic file
  `build/java-reference/source-vandy-atlas-first-loop.log`.
- Python log after changing the production atlas loader: ignored diagnostic
  file `build/java-reference/python-source-vandy-postfix-first-loop.log`.
- Focused regression command after the implementation change:
  `PYTHONPATH=src .venv/bin/pytest -q tests/test_atlas_assets.py tests/test_toads_atlas_startup.py tests/test_toads_registered_first_iteration_cycle.py`
  — **14 passed in 9.64s**.
- The deterministic inputs are committed under
  `tests/data/synthetic/toads_downsampled_vandy_atlas/`; the atlas NIfTIs are
  bundled under `data/atlases/`.
- The low-memory 12-iteration Python diagnostic completed in 37.7 seconds
  with history and per-event capture disabled. A trace retaining all thinning
  and growing events exited with code 137; the no-history run yielded the
  same algorithm settings and final classification.
- No full pytest suite was run for this change.

## Twelve-iteration compact comparison

Both Java and Python completed 12 iterations on this compact synthetic case,
but neither reached the `1e-4` energy-ratio threshold. Java ended at
`0.02333523`; Python ended at `0.0003040527`. Their ratios track closely
through iteration 2, then diverge; by iteration 12 the objective values and
ratios no longer match. This keeps convergence/source-equivalence work open.

Despite that internal difference, their final native-grid segmentations
agree on 54,981 of 55,080 voxels (99.82%). Every one of the 11 labels has Dice
at least `0.9827` (the minimum is label 22); all meet the requested 0.97
threshold for this compact TOADS-only synthetic result. This is not the full
MaCRUISE output acceptance result because Dura, CRUISE, and final MaSeg
refinement are not part of this pairwise comparison.

The final Python NIfTI is `(36,45,34)` uint8 with RAS axes. Its terminal
classification hash is recorded in the ignored diagnostic log
`build/java-reference/python-source-vandy-12iter-low-memory.log`; Java's
12-iteration completion is in
`build/java-reference/source-vandy-atlas-12iter.log`. These arrays/logs are
build artifacts and are not committed. No full pytest suite was run for this
change.

## Brain-resolution source-order TOADS run

After the source-order loader change, the bounded `toads-loop` CLI completed
against the generated `(172,216,164)` synthetic T1/mask and the bundled
audited Vandy atlas descriptor. To regenerate the ignored large synthetic
input directory and repeat the run from the repository root:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise toads-loop \
  --t1-n4 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output build/vandy-brain-resolution-memory/toads-source-order-two-iterations.nii.gz \
  --maximum-iterations 2 --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --padding 10 --relation-limit 3 \
  --provenance build/vandy-brain-resolution-memory/toads-source-order-two-iterations.json
```

The command exited 0 in **525.4 seconds**. The RSS watchdog sampled a peak of
**6.232 GB** and did not terminate it. TOADS stopped at the requested
two-iteration limit with terminal ratio `0.323834717`; convergence was not
reached. The output is finite uint8 NIfTI on the `(172,216,164)` RAS grid with
labels `[1,5,6,14,15,16,17,18,22,24,25]`. Output SHA-256:
`736e1796276a09ffbfee7ade0affe8fb7e2d359a28e83bedc0f5efb5a649fc4f`.
The JSON provenance records the bundled atlas descriptor and effective
parameters. The `.nii.gz`, JSON, and full log are generated in ignored
`build/vandy-brain-resolution-memory/`; this is a TOADS-only synthetic run,
not a completed base workflow or real-data validation.

Next, localize the post-iteration-2 divergence and rerun the restartable base
workflow using source-ordered atlas inputs. Keep these compact synthetic
results separate from brain-resolution convergence and whole-workflow
acceptance.

## Brain-resolution source-order base workflow

The source-order atlas path has now completed the six-stage synthetic base
workflow on generated `(172,216,164)` brain-resolution T1, brain-mask, and
MaSeg inputs. The run used the bundled, audited Vandy descriptor and bounded
stage settings (two TOADS iterations per pass, one Dura iteration, one ACE
diffusion iteration, and one TGDM iteration per surface). This validates stage
composition, output-grid checks, manifest hashing, and restart behavior. It
does not establish TOADS convergence, a whole-workflow Java comparison, Dice
acceptance, or real-data validation.

Reproduce from the repository root:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --maseg build/vandy-brain-resolution-memory/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/vandy-brain-resolution-memory/source-order-base-run-reproduction \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

The workflow exited 0 under the 16 GB RSS watchdog. It completed all six
stages and their manifest checks: N4, first TOADS, Dura, second TOADS, CRUISE
core, and MaSeg refinement. Peak sampled RSS was **6.427 GB**. A second run
with the same command plus `--resume` verified and skipped all six stages,
exited 0, and peaked at **0.126 GB**. The complete manifest and products are
under ignored `build/vandy-brain-resolution-memory/source-order-base-run/`.

The manifest is `status=complete`, `converged=false`, and labels the workflow
scope as `synthetic/native-grid base route; not a full Java-equivalence claim`.
Both TOADS stages stopped at their two-iteration caps without reaching the
energy threshold: the first terminal ratio was **1.10402799**, and the second
was **0.571384311**. All three GVF components emitted a warning that they did
not converge in 200 iterations. These are explicit outstanding algorithmic
limits, despite successful stage completion.

Every stage output passed its grid/finite-value/hash check. The final inner,
central, and outer level sets are finite RAS NIfTIs on the original
`(172,216,164)` grid; `SEG_refine.nii.gz` is also finite and RAS on that grid.
The CRUISE internal crop is `(163,207,156)` RAS, and the 3-component GVF is
`(163,207,156,3)` RAS. The manifest SHA-256 values for the final level sets
are, respectively, `9ede89a4...e968a2` (inner), `ba537852...bd92f4c`
(central), and `694d65c6...2b1c378` (outer); refined segmentation is
`112d2d2a...118dc124`. Full hashes and stage output records are retained in
the run manifest. The input fixture, stage products, and watchdog logs are
ignored build artifacts, not committed test samples.

The production loader behavior and focused regression evidence are on code
commit `4b57532`: `PYTHONPATH=src .venv/bin/pytest -q
tests/test_atlas_assets.py tests/test_toads_atlas_startup.py
tests/test_toads_registered_first_iteration_cycle.py` — **14 passed in
9.64s**. No full pytest suite was run for this integration. The most recent
full-suite evidence remains **1,138 passed, 0 failed, 0 skipped, 11 expected
RAS warnings in 793.64s** on `ad34321`, before the source-order loader
change.

## GVF convergence warning and selected-source probe

The full-resolution CRUISE-core run produced one SciPy CG iteration-cap warning
for each vector component at 200 iterations and relative tolerance `1e-6`.
Applying the same selected source system to the stored full-resolution output
gave relative residuals `3.2041e-4`, `3.1943e-4`, and `3.5092e-4` for x/y/z.
The warning is therefore genuine for the Python stopping rule. It should not
be described as failure of Java convergence: the selected VandyCRUISE
`Poisson3d` runs the requested fixed number of multigrid V-cycles and does not
use this CG residual stopping test. That solver substitution has now been
compared directly on the full `(163,207,156)` crop; it fails the existing
numerical gate substantially at brain resolution.

To measure the solver substitution on the same brain-resolution-derived ACE
input with the same 200-iteration setting, I ran a stride-4 comparison against
the copied Java harness. The probe is `(41,52,39)`, uses lambda 1.0, and
compares each output component. All values are finite; each component meets
the existing GVF synthetic gate of maximum absolute error ≤`5e-5`:

| Component | MAE | Maximum absolute error | RMSE |
|---|---:|---:|---:|
| x | `6.8446e-6` | `3.2157e-5` | `9.7091e-6` |
| y | `2.6879e-6` | `2.2449e-5` | `4.2178e-6` |
| z | `2.6789e-6` | `3.6813e-5` | `4.8954e-6` |

Reproduce the comparison after running the base workflow:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/compare_gvf_java_probe.py \
  --input build/vandy-brain-resolution-memory/source-order-base-run/stages/cruise-core/attempt-01/ace_enhanced_gm.nii.gz \
  --factor 4 --iterations 200 --lambda-value 1.0 \
  --output-dir build/java-reference/vandy-gvf-source-order-probe-reproduction
```

The script compiles the selected upstream GVF sources and stores the
comparison JSON plus temporary Java binary interchange files in the ignored
build directory. This probe checks the solver substitution on a subsampled
input; it does not validate the full-resolution field numerically. Retain the
full-resolution warning and residual in workflow provenance until a direct
Java reference comparison or a more faithful multigrid port resolves that
gap. No full pytest suite was run for this probe.

### Full-grid Java comparison

The same probe was run at `--factor 1` with the complete `(163,207,156)` ACE
enhanced GM volume, lambda 1.0, and 200 iterations for each implementation.
The Java multigrid reference completed in about one minute. The combined
watchdog sampled **1.053 GB** peak RSS against the 16 GB cap; both solvers
returned finite vectors. Their full-grid numerical differences are far above
the existing synthetic gate:

| Component | MAE | Maximum absolute error | RMSE |
|---|---:|---:|---:|
| x | `2.1823e-3` | `2.0948e-2` | `3.6712e-3` |
| y | `2.4411e-3` | `2.1780e-2` | `4.0976e-3` |
| z | `2.7043e-3` | `2.1898e-2` | `4.4090e-3` |

Reproduce with the same command and `--factor 1`, writing to a distinct
ignored build directory. Full output and solver warnings are in
`build/java-reference/vandy-gvf-source-order-full/comparison.json` from the
recorded run. The stride-4 agreement is not representative of the full-grid
system. The CG implementation cannot be considered a safe source substitute
for brain-resolution GVF; the next implementation step is to port the
selected VandyCRUISE multigrid V-cycle, then compare it against the Java
harness at both compact and full crop resolution.

## Source multigrid port and full-grid verification

`generalized_gvf` now uses a Numba-compiled port of the selected Java
`Poisson3d` V-cycle instead of SciPy CG. It preserves the Java coefficient
and residual 2x2x2 downsampling, nearest-neighbor transfer, coarse-grid
correction, clamped-boundary residual, spacing factors, and eight-octant
red/black Gauss-Seidel update order. The default stopping criterion is the
source's absolute max residual `1e-6`, with the selected 200-cycle maximum.
The `solver_tolerance` API/CLI parameter now names this absolute residual
threshold; the source default is unchanged.

The stride-4 brain-derived comparison reproduces with the same probe script
using the updated Python implementation. Its maximum component errors are
`5.07e-7`, `4.92e-7`, and `6.41e-7` for x/y/z.

The full-grid direct comparison was repeated on the actual `(163,207,156)`
RAS ACE-enhanced GM crop, at lambda 1.0 and 200 iterations/cycles:

| Component | MAE | Maximum absolute error | RMSE |
|---|---:|---:|---:|
| x | `9.6787e-8` | `6.7651e-6` | `2.2850e-7` |
| y | `4.4645e-8` | `4.3958e-7` | `7.2528e-8` |
| z | `4.3787e-8` | `5.8115e-7` | `7.2254e-8` |

All three components pass the existing max-error gate `5e-5`, and both
outputs are finite. The 16 GB watchdog sampled a combined peak of **0.737
GB** and the comparison exited 0. Reproduce after the generated base fixture
exists with:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  env PYTHONPATH=src ./.venv/bin/python scripts/compare_gvf_java_probe.py \
  --input build/vandy-brain-resolution-memory/source-order-base-run/stages/cruise-core/attempt-01/ace_enhanced_gm.nii.gz \
  --factor 1 --iterations 200 --lambda-value 1.0 \
  --output-dir build/java-reference/vandy-gvf-multigrid-reproduction
```

The runnable report JSON from this validation is under ignored
`build/java-reference/vandy-gvf-multigrid-full/comparison.json`. This is a
full-grid, selected-component Java comparison; it still does not establish
full MaCRUISE source equivalence.

### Source-multigrid CRUISE-core and base-run resume

The new solver was then run through `toads-cruise-core` on the same full
resolution second-pass TOADS products. It produced all eleven scientific
outputs plus provenance in RAS; its 16 GB watchdog sampled a **2.304 GB**
peak. The final inner, central, and outer level sets remained very close to
those from the prior CG-based stage: sign-mask Dice was 1.0, 0.999999, and
0.999986, respectively; symmetric mean nearest-surface distances were
`3.29e-6`, `9.44e-6`, and `2.95e-4` mm; Hausdorff distances were 0.077, 0.416,
and 0.485 mm. These compare two Python stages, not either result to Java
MaCRUISE surfaces.

MaSeg refinement from the new inner/outer level sets changed 17 voxels from
the previous CG-based refined result. Across labels `[0,1,4,100]`, minimum
Dice was **0.9999436**; the output is uint8, `(172,216,164)`, and RAS.

Finally, the base workflow's `run.json` was preserved as
`run.before-source-multigrid.json`, the old core stage was marked stale, and
the restartable base command was resumed with unchanged inputs/settings. It
verified/skipped N4, first TOADS, Dura, and second TOADS; recomputed core as
attempt 2 and refinement as attempt 2; and passed the stage checks for all 13
core outputs and all 3 refinement outputs. The manifest is again
`status=complete`, `converged=false`. TOADS still hit its two-iteration caps
at ratios `1.10402799` and `0.571384311`; this integration validates the
source-multigrid core and restart path, not TOADS convergence or whole-workflow
Java equivalence. The resume run peaked at **2.474 GB** under the 16 GB cap
and exited 0. Final source-multigrid RAS outputs are in the manifest's
`cruise-core/attempt-02` and `maseg-refinement/attempt-02` records.

Focused regression command after replacing CG:
`PYTHONPATH=src ./.venv/bin/pytest -q tests/test_gvf.py tests/test_plus_gvf.py`
— **8 passed in 1.49s**. Full-suite evaluation is pending. The earlier CG
mismatch above is preserved as the diagnostic that motivated the change; the
source-multigrid results supersede its statement that the solver substitution
remains a blocker.
