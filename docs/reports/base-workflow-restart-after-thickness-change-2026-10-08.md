# Synthetic base workflow restart after thickness changes — 2026-10-08

## Reproduction

The six-stage synthetic base workflow was rerun on the committed `(172,216,164)` RAS fixture with a 16 GB sampled-RSS watchdog:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/ace-indexed-heap-2026-10-07/base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1 --resume
```

This command rebuilt invalidated stages in a new attempt generation while preserving prior attempts. It completed with exit code 0, all six stage checks passed, the manifest reports `status=complete` and `converged=true`, and 16 final products are recorded. The run contains 38 hashed stage outputs across the six stages; its NIfTI outputs are finite and RAS. First TOADS converged at iteration 12 with terminal energy ratio `6.96809075e-05`. Second TOADS converged at iteration 34 with terminal ratio `5.54052858e-05`. Both traces and terminal ratios match the previously recorded selected-Java comparison. Peak sampled RSS was **8.190 GB**; the 16 GB watchdog did not terminate the child.

The same command and settings were immediately rerun with `--limit-gb 4`. Resume verified and skipped N4, first TOADS, Dura, second TOADS, CRUISE core, and MaSeg refinement; exit code was 0 and peak sampled RSS was **0.123 GB**. This verifies the current manifest and restart behavior.

## Focused test evidence

On implementation/documentation commit `b500127`, before this evidence-only report update:

```sh
./.venv/bin/python -m pytest -q tests/test_base_workflow_restartable.py
```

Result: **3 passed, 0 failed, 0 skipped in 41.40s**. No full pytest suite was run in this cycle. This is synthetic integration and selected-stage source-comparison evidence, not whole-workflow MATLAB/JIST equivalence or real-MRI validation.

The base CLI writes cropped tissue and ACE products, GVF, inner/central/outer level sets, and refined MaSeg labels. It does not create meshes or thickness maps. The selected JIST thickness wrapper comparison remains a separate held-fixed-mesh source comparison; this workflow rerun does not validate mesh extraction or thickness integration.

The supplied Singularity output directory remains unavailable as a comparison source, and the available human T1 still requires a matching brain mask. Whole-workflow source equivalence, refinement comparison, mesh-connectivity differences, plus-workflow parity, and real-MRI validation remain open.
