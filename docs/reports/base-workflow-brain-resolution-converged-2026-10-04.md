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

The final inner, central, outer, and refined MaSeg NIfTI outputs are finite and
on the native `(172,216,164)` grid with RAS axis codes. The run's peak sampled
RSS was **8.308 GB** under a **16 GB** process watchdog; it exited with code 0
and was not terminated. A same-parameter `--resume` invocation verified and
skipped all six stages, exited with code 0, and peaked at **0.122 GB**.

## Reproduction

From the repository root, regenerate the deterministic atlas-derived phantom
in an ignored build directory:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
```

Run the complete synthetic base workflow with converged TOADS iteration limits:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --maseg build/vandy-brain-resolution-memory/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/vandy-brain-resolution-memory/source-order-base-run-converged \
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
