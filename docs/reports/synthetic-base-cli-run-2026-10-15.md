# Restartable synthetic BASE CLI run — 2026-10-15

## Runnable command

From the repository root, the documented `base-run` command executes the
compact synthetic T1, supplied brain mask, bundled atlas descriptor, and
MaSeg inputs. Its MaSeg must contain some cortical gray matter because the
selected BASE crop is defined from the MaSeg-derived GM map. The fixture
generator now includes 125 label-100 voxels as crop support, in addition to
the eight brainstem/cerebellum cleanup labels.

```bash
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/test-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 5 --report-interval-seconds 60 -- \
  .venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
  --output-dir build/restartable-base-command-2026-10-15-rerun \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 2 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

The command completed at **0.599 GB peak sampled RSS**. All eight manifest
stages passed validation: N4, first TOADS, Dura, second TOADS, base tissue
preparation, CRUISE core, MaSeg refinement, and surface/thickness. The manifest
reports 24 final outputs. Repeating the command with `--resume` verified and
skipped all eight stages; no stage reran.

Both TOADS passes stop at the explicit two-iteration cap, with energy ratios
`0.0556259789` and `0.118235528`; neither met its convergence threshold. This
is a completed restartable synthetic execution, not a convergence or
source-equivalence claim. The focused restartable workflow tests passed **5
tests, 0 failed, 3 expected orientation warnings in 14.54 seconds** under an
8 GB watchdog; peak sampled RSS was **0.639 GB**. The full suite was not run.
