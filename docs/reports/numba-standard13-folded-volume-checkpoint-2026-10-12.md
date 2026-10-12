# Ordered standard-13 acceleration and reproducible PLUS checkpoint

## Source-behavior implementation

The standard-13 TOADS membership sweep now has a Numba implementation and retains the Python implementation as an executable oracle (`backend="python"`). Numba is the default. The accelerated kernel preserves source traversal and in-place neighbor updates, class and energy state ordering, float32 arithmetic, the six-face neighborhood, and source-selected zero-variance behavior.

The accelerated and Python backends match bit-for-bit on the 5³ Java-harness case and on a 15³ folded synthetic volume tiled from that input. All membership, relation, energy, denominator, objective-history arrays, and the scalar objective agree exactly. The focused group passed **11 tests in 6.49 s** on implementation commit `e09e104c`. This establishes Python-backend equivalence on synthetic states; it is not a new Java-equivalence claim.

## ANTs reproducibility correction

Before the correction, repeated Rigid registrations on identical synthetic inputs and settings produced different transforms. With the PLUS sampling rate 0.25, two runs differed at **100,955 voxels** (maximum absolute difference 0.41469). ANTsX's reproducibility guidance calls for a fixed random seed and single-threaded execution; setting those variables only after ANTsPy/ITK had initialized was too late. PyMaCRUISE now sets `ANTS_RANDOM_SEED` and `ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1` before converting images to ANTs objects, then restores the prior environment.

After the fix, repeated registrations produced bit-identical transforms and warped images. The committed regression test repeats Rigid registration with fixed seed and verifies exact equality. Registration and CLI tests passed **9 tests in 4.88 s** on commit `bca4e8cc`. The source guidance is [ANTsX's reproducibility notes](https://github.com/ANTsX/ANTs/wiki/antsRegistration-reproducibility-issues) and [ANTsPy's environment-variable notes](https://github.com/ANTsX/ANTsPy/wiki/Important-environment-variables).

## Restartable synthetic PLUS run

Runnable command (16 GB RSS watchdog; single-threaded numerical libraries):

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 30 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_labeled_folded_cortex/t1.nii.gz \
  --flair tests/data/synthetic/plus_labeled_folded_cortex/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_labeled_folded_cortex/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_labeled_folded_cortex/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-labeled-folded-deterministic-numba-onepass-2026-10-12 \
  --max-passes 1
```

All **nine workflow stages completed** under the watchdog; the manifest is `build/plus-labeled-folded-deterministic-numba-onepass-2026-10-12/plus-full-workflow.json`. It is a one-pass checkpoint with `terminal=false`, not a converged result. The iteration-2 membership update took **1.065 s**; registered shape-prior computation took 222.47 s and is now the larger measured bottleneck. A prior unseeded 65³ run recorded about 1,014 s for the membership stage, but that is not a controlled speed comparison because its registration-derived inputs differed.

The current run's registered FLAIR output is bit-identical to two independently repeated registrations with the same seed and settings. Its first ribbon audit still fails: GM-between-surfaces support **82.98%**, WM-inside-inner support **89.92%**, and closest-surface thickness p05/median/p95 **0.174/1.887/5.615 mm**. The level-set mean thickness is 4.690 mm, but nested voxel regions fail. All meshes are watertight but each has 21 connected components and Euler characteristic 36. These are synthetic checkpoint measurements, not Java equivalence or real-human PLUS validation.

## Next integration work

Keep the Python backend as an oracle and the deterministic ANTs setup enabled. Diagnose the synthetic ribbon support, nesting, and mesh-component failures, then resume the checkpoint for additional source-selected passes and repeat the explicit quality audit. A held-fixed full-workflow Python-backend comparison has not yet been run. The previous whole-run comparison is no longer used to assess backend equivalence because it predates repeatable registration.
