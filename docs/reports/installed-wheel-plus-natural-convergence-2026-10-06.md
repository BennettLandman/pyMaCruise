# Installed-wheel plus natural-convergence run — 2026-10-06

## Result

Built PyMaCRUISE from the committed checkout and installed the wheel into an
isolated target. The process ran from
`build/wheel-plus-natural/`, outside the source checkout, and imported
`pymacruise` from `build/wheel-plus-natural/site/pymacruise`. It used the
committed 36×45×34 folded RAS T1/FLAIR/user-mask/MaSeg fixture and the
committed compact synthetic LesionTOADS atlas.

All nine stages completed. The top-level manifest reports `complete` and
`terminal=true`; the selected defaults were tolerance `0.0001`, 99 maximum
passes, and 99 maximum iterations. LesionTOADS stopped by convergence at
iteration 5 with energy ratio `8.241410796472337e-06`. All nine child-stage
manifests and their 23 recorded outputs passed SHA-256 verification; every
NIfTI output was finite and RAS. The final-output hashes also verified. Surface
thickness reports a mean of `0.11254296722799782 mm` over 998 vertices, and
the volume thickness product is available with one crossing voxel. A second
invocation with the same options verified and reused all nine stages.

Peak sampled RSS for the initial run was **10.229 GB** under the 16 GB
watchdog. The span between the first and last stage-manifest modification
times was **3,229.14 seconds (53 minutes 49 seconds)**. The resumed command
completed with peak sampled RSS **0.126 GB**. Fontconfig printed warnings
about unavailable user cache directories in this restricted environment;
the workflow completed and all output checks passed.

This is installed-wheel synthetic execution and convergence evidence. It is
not Java/MATLAB whole-workflow equivalence, source numerical acceptance,
brain-resolution plus validation, or real-MRI validation.

## Comparison with the source-checkout run

I compared the installed-wheel outputs with the P41v source-checkout test on
the same committed synthetic inputs. The checkout test used an iteration
ceiling of 50; the wheel used the default ceiling of 99. Both stopped at
iteration 5 with the same energy ratio. Discrete LesionTOADS segmentation,
classification, and subtype arrays matched exactly. All three CRUISE level
sets, the volume-thickness array, the extracted surface meshes, and the
thickness summary also matched exactly.

Some earlier float products were not byte-identical across these separate
runs: maximum absolute differences were `1.52587890625e-05` in the
preprocessed modality volume, `1.7881393432617188e-07` in recropped
modalities, and `1.9073486328125e-06` in initialized LesionTOADS memberships.
The cause of this small upstream variation has not been isolated. This
comparison therefore supports exact final synthetic segmentation/surface
agreement for these two runs, but does not establish bitwise equality of all
intermediate files or source-algorithm equivalence.

## Reproduction

Build and install (the wheel was built from `981e708`; the implementation is
unchanged from the full-suite-tested commit `09a7bd0`):

```bash
mkdir -p build/wheel-plus-natural/dist build/wheel-plus-natural/site
./.venv/bin/python -m build --wheel --no-isolation \
  --outdir build/wheel-plus-natural/dist .
./.venv/bin/python -m pip install --no-deps --target \
  build/wheel-plus-natural/site \
  build/wheel-plus-natural/dist/pymacruise-0.1.0-py3-none-any.whl
```

Run from `build/wheel-plus-natural/`, with `PYTHONPATH` pointing only to the
isolated wheel target and `PYTHONNOUSERSITE=1`. This command uses the committed
fixture files and `plus-full-run`; registration/N4/first-loop cost settings
are reduced for synthetic runtime, while convergence tolerance, pass ceiling,
and iteration ceiling remain at their CLI defaults:

```bash
env PYTHONNOUSERSITE=1 PYTHONPATH="$PWD/site" \
  MPLCONFIGDIR="$PWD/mpl-cache" \
  /absolute/path/to/pymacruise/.venv/bin/python \
  /absolute/path/to/pymacruise/scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  /absolute/path/to/pymacruise/.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 /absolute/path/to/pymacruise/tests/data/synthetic/plus_larger_grid/t1.nii.gz \
  --flair /absolute/path/to/pymacruise/tests/data/synthetic/plus_larger_grid/flair.nii.gz \
  --brain-mask /absolute/path/to/pymacruise/tests/data/synthetic/plus_larger_grid/brain-mask.nii.gz \
  --maseg /absolute/path/to/pymacruise/tests/data/synthetic/plus_larger_grid/maseg.nii.gz \
  --atlas-descriptor /absolute/path/to/pymacruise/tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir natural-run \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

Append `--resume` to that same command for the reuse check. Both the initial
and resume commands ran under the watchdog with a 16 GB limit.

The full suite on the unchanged implementation commit
`09a7bd034f1d32c336525a9d50002af897f79140` passed **1,219 tests, 0 failed,
0 skipped, 11 expected RAS warnings in 1,699.06 seconds** with peak sampled
RSS **13.052 GB**; see [the full-suite report](plus-natural-convergence-2026-10-06.md).
