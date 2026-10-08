# Current-code synthetic base workflow and restart — 2026-10-08

## Scope

This run exercises the full seven-stage synthetic base workflow on the brain-resolution fixture after the VandyCRUISE rigid-pyramid source correction (`069e9d3`, implementation signature `abe8c0ecb3a5127e7540d515d345ca425f458d35508f00eeab4e9b722611f457`). It is restartable synthetic integration evidence, not whole-workflow MATLAB/JIST equivalence or real-MRI validation.

## Run and result

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.2 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/jist18-base-workflow-2026-10-08/run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1 \
  --resume
```

The pre-existing manifest belonged to an older implementation signature. With `--resume`, the runner preserved old attempts and regenerated stages affected by the new signature. Exit code was 0 and all seven stage checks passed. Both TOADS loops reached their energy thresholds (not iteration caps): first pass, 12 iterations and terminal ratio `6.96809075e-05`; second pass, 34 iterations and terminal ratio `5.54052858e-05`. The second pass also produced 544,793 changed classification voxels, matching the previous successful run. The run manifest reports all seven stages complete.

An independent audit verified all **46 stage-output SHA-256 hashes**, loaded all **31 stage NIfTIs** to confirm finite voxel values and RAS orientation, and found the four expected mesh outputs. The 16 GB RSS watchdog sampled a **7.690 GB peak** and did not terminate the run.

## Restart check

The same command and output directory with `--resume` exited 0 and printed `Resume: verified ...; skipping` for each of the seven stages (N4, TOADS-first, Dura, TOADS-second, CRUISE core, MaSeg refinement, and surface-thickness). The resume watchdog sampled **0.124 GB peak RSS**. No stage was recomputed.

## Interpretation and limits

This confirms current-code stage composition, convergence metadata, output integrity/orientation, and restart behavior for this deterministic synthetic case. It does not establish full MaCRUISE source equivalence. Selected Java comparison for the updated VandyCRUISE pyramid is documented separately in [the actual-source registration report](toads-registration-pyramid-actual-source-2026-10-08.md); the direct method harness and focused tests there are not substitutes for full-workflow comparison. Real MRI remains unvalidated because the human T1 needs its matching brain mask.
