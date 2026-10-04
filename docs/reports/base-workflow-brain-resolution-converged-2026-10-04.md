# Brain-resolution synthetic base workflow integration — 2026-10-04

## Result

The restartable six-stage `base-run` completed on a generated brain-resolution
synthetic input set `(172,216,164)` using the bundled VandyCRUISE atlas. N4,
both TOADS passes, Dura, CRUISE core, and MaSeg refinement all passed their
stage-output checks. The first TOADS pass reached the 50-iteration cap without
meeting its energy threshold (terminal ratio `0.0301428866`); the second pass
converged at iteration 20 (ratio `9.09511073e-05`). Consequently, the workflow
manifest says `status=complete` but `converged=false` and
`source_equivalence_claim=false`. This is successful end-to-end synthetic
integration evidence, not a converged full workflow, Java equivalence, or
clinical validation.

## Source-configured stopping rule

The MATLAB V3 wrapper `run_Toads_jist.m` configures the selected TOADS loop
with `-inMaximum 1.0E-4` and `-inMaximum2 50` (source SHA-256
`8c2e9e8ea554eaf4d731082ae1c6d30e66130d6ccf048e5373d622783eb52575`). Thus
the integrated first pass stopping at iteration 50 with its energy threshold
false is consistent with the wrapper's configured iteration cap. It is not,
on its own, evidence of a Python/source mismatch. The second pass stops earlier
because it reaches the configured threshold. The manifest retains
`converged=false` because the first pass did not meet the threshold; stage
completion means the configured loop returned and the workflow continued, not
that both loops satisfied the energy threshold. The remaining question is
whether the selected Python and Java implementations produce sufficiently
similar intermediate and final results under these same source settings.

## Standalone versus integrated first-pass input

The standalone brain-resolution convergence report used the generated raw T1
fixture as the `--t1-n4` input. The base workflow's first pass uses the N4
corrected output instead. The atlas descriptor used in the integrated run is
byte-identical to the bundled source descriptor (SHA-256
`cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939`), and the
registration iteration controls and 50-iteration cap match the standalone
command. The raw and N4 images are on the same grid, but have different
intensity distributions within the supplied brain mask:

| Input | Mean | Standard deviation | Minimum | Maximum |
|---|---:|---:|---:|---:|
| Generated raw T1 | 444.8967 | 230.8274 | 152.48 | 788.00 |
| N4-corrected T1 | 366.6892 | 189.6747 | 113.4732 | 689.9545 |

The standalone and integrated energy trajectories are therefore not a
controlled equivalence comparison. In the integrated first pass, iteration 1
was `0.161843687`, iteration 2 was `1.10402799`, and iteration 50 ended at
`0.0301428866`; the standalone raw-T1 run had iteration-1 ratio `1.66407061`
and converged at iteration 20. This identifies the input stage as a condition
to isolate, but does not establish that N4 alone caused the different
trajectory or that either trajectory matches Java/MATLAB. The integrated
first-pass output is direct evidence for the N4-input route; it should be
compared with a source-matched workflow reference before changing N4 or TOADS
behavior.

An additional raw-versus-N4 experiment used the committed compact
`toads_downsampled_vandy_atlas` phantom `(36,45,34)`, its generated atlas
descriptor, one shared brain mask, and otherwise matching 50-iteration,
padding-10 TOADS settings. The raw-T1 run reached its cap at ratio
`0.00278222049`; the N4-input run converged at iteration 20 at
`6.68482171e-05`. Their classifications differed in 1,488 of 55,080 voxels
(97.30% voxel agreement; minimum per-label Dice `0.748344`, label 17). This
compact A/B is synthetic diagnostic evidence only. Its direction does not
explain the opposite brain-resolution result, showing that this response is
fixture/scale-dependent; it does not prove N4 causes or fixes convergence.
The voxel agreement and Dice compare the raw run stopped at iteration 50 with
the N4 run stopped at iteration 20; these terminal-output differences are
descriptive and are not a numerical acceptance result.

The compact A/B can be reproduced from the committed
`tests/data/synthetic/toads_downsampled_vandy_atlas/` fixture. First create
the N4 input:

```bash
./.venv/bin/pymacruise n4 \
  --input tests/data/synthetic/toads_downsampled_vandy_atlas/t1.nii.gz \
  --output build/toads-n4-convergence-diagnostic/t1-n4.nii.gz
```

Run each arm with the same 50-iteration settings (change `--t1-n4` and the
output/provenance names for the N4 arm):

```bash
./.venv/bin/pymacruise toads-loop \
  --t1-n4 tests/data/synthetic/toads_downsampled_vandy_atlas/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_downsampled_vandy_atlas/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_downsampled_vandy_atlas/atlas.txt \
  --output build/toads-n4-convergence-diagnostic/raw-classification.nii.gz \
  --maximum-iterations 50 --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --padding 10 --relation-limit 3 --maximum-distance 1e-4 \
  --provenance build/toads-n4-convergence-diagnostic/raw.json
```

The final inner, central, outer, and refined MaSeg NIfTI outputs are finite and
on the native `(172,216,164)` grid with RAS axis codes. The run's peak sampled
RSS was **8.308 GB** under a **16 GB** process watchdog; it exited with code 0
and was not terminated. A same-parameter `--resume` invocation verified and
skipped all six stages, exited with code 0, and peaked at **0.122 GB**.

## Reproduction

From the repository root, regenerate a comparison copy of the committed
brain-resolution inputs in an ignored build directory:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
```

Run the complete synthetic base workflow with the wrapper's configured TOADS
stopping rule (energy threshold `1e-4`, iteration cap 50):

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/vandy-brain-resolution-memory/committed-input-base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Append `--resume` to that command to verify the cached outputs and skip completed
stages. This generated fixture is deliberately synthetic and is not a substitute
for a Java/MATLAB whole-workflow reference.

## Output evidence

`run.json` records all six stages as complete, output hashes, RAS axis codes,
the two TOADS stop records, `converged=false`, and
`source_equivalence_claim=false`. Representative final-output SHA-256 values:

| Output | Shape | Axis codes | SHA-256 |
|---|---:|---|---|
| `inner.nii.gz` | `(172,216,164)` | RAS | `b2c3f6fee3435fb9d6dd08c698214f062d09979f1b82615e41c151b2da6930e8` |
| `central.nii.gz` | `(172,216,164)` | RAS | `cb282838f14c127dc0bf70ba0630b9b3cba299c13bef062ed360eecd6aa69661` |
| `outer.nii.gz` | `(172,216,164)` | RAS | `985c228de883de14328fe634ccdf0f12c8c9eff014e2723785b338bea246e0e7` |
| `SEG_refine.nii.gz` | `(172,216,164)` | RAS | `c3f6aed851880fe2bc9a3de95df2d1a04b182325bd2117bf226d5a599ec8b27c` |

No pytest suite was run for this runtime experiment. The earlier full suite on
`ea31121` remains the latest full-suite evidence and does not cover this
integration run.
