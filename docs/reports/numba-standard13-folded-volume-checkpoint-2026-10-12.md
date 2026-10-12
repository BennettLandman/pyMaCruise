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

### Two-pass continuation and restart recovery

The one-pass checkpoint was resumed after an interrupted attempt. The runner recorded the incomplete attempt as interrupted, reused the verified preprocessing/crop/alignment/registered-branch stages, and completed the nine-stage workflow. This verifies restart recovery for that interruption boundary. The resulting manifest is `build/plus-labeled-folded-deterministic-numba-twopass-2026-10-12/plus-full-workflow.json`; all nine public stages are complete, but the lesion loop remains nonterminal at iteration 3. The run peaked at **7.124 GB sampled RSS**, below its 16 GB watchdog.

The controlled comparison holds the registered upstream state fixed: registered FLAIR is bit-identical, and all **79,613,040** registered-branch membership values are bit-identical between the one- and two-pass runs. Continuing the lesion loop changes **10,826,827 / 79,613,040** registered-loop membership values (mean absolute change 0.005218; maximum 1.0) and **119,918 / 4,683,120** segmentation voxels. The later surface audit is worse on its principal label-support measure: GM support falls from **82.98%** to **73.59%**; WM support is **89.81%** and CSF support remains 100%. Two-pass meshes have 24 connected components each (watertight, Euler 32); nesting fails by one inner-outside-central voxel. Production closest-surface thickness p05/median/p95 is **0.160/1.395/5.582 mm**, and level-set mean thickness is 4.545 mm. The latter passes its broad mean bound, but support, topology, nesting, and production thickness checks fail. The loop's energy ratio is NaN under the selected source zero-variance behavior and it has no terminal stop reason.

This is a negative synthetic result: an additional continuation pass does not repair the ribbon and makes GM support worse. Since registration and the registered-branch membership input match exactly, the quality change is downstream in the continuing lesion loop and its surface products. It does not by itself identify a Python defect: continued state evolution is expected, and there is no held-fixed Java replay for this run. These results remain synthetic workflow evidence only, not source equivalence or in-vivo PLUS validation.

## Next integration work

Keep the Python backend as an oracle and the deterministic ANTs setup enabled. Compare the one-pass and two-pass continuation states against selected Java methods at the same saved state, then isolate the first state transition associated with support loss. Do not spend additional full-run cycles merely increasing the pass count until this diagnosis identifies a meaningful experiment. A held-fixed full-workflow Python-backend comparison has not yet been run. The previous whole-run comparison is no longer used to assess backend equivalence because it predates repeatable registration.
