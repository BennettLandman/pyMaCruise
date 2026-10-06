# Installed-wheel synthetic MaCRUISEplus run — 2026-10-05

## Result

Built and installed the PyMaCRUISE wheel, then ran the installed package's
`plus-full-run` CLI from a working directory outside the source checkout. The
committed 25³ folded-cortex fixture recipe in
[`tests/test_plus_full_workflow.py`](../../tests/test_plus_full_workflow.py)
created T1, shifted-grid FLAIR, the matching brain mask, and MaSeg inputs. The
run used the committed compact synthetic LesionTOADS descriptor and its NIfTI
assets as an explicit input.

The run completed all nine stages:

1. modality registration and N4 preprocessing
2. crop
3. initial atlas alignment
4. registered recrop
5. terminal LesionTOADS loop
6. Dura outputs
7. native-T1 transfer
8. MaSeg cortical separation
9. CRUISE surface and thickness workflow

`plus-full-workflow.json` reports `status=complete` and `terminal=true`. Every
stage manifest and output hash verified; all NIfTI outputs were finite and
RAS. The surface thickness summary reports a finite mean of `0.06391155 mm`
over this synthetic case. A second invocation using `--resume` verified and
reused the completed run. Peak sampled RSS was **1.801 GB** under the 16 GB
watchdog.

## Reproduction

Build and install the wheel in an isolated target directory:

```bash
./.venv/bin/python -m build --wheel --no-isolation \
  --outdir build/wheel-plus-run/dist .
mkdir -p build/wheel-plus-run/site
./.venv/bin/python -m pip install --no-deps --target \
  build/wheel-plus-run/site \
  build/wheel-plus-run/dist/pymacruise-0.1.0-py3-none-any.whl
```

From `build/wheel-plus-run`, run the installed package with inputs generated
by [`_inputs`](../../tests/test_plus_full_workflow.py) and options matching
[`_command`](../../tests/test_plus_full_workflow.py):

```bash
PYTHONNOUSERSITE=1 PYTHONPATH="/absolute/path/to/pymacruise/build/wheel-plus-run/site" \
  MPLCONFIGDIR="/absolute/path/to/pymacruise/build/wheel-plus-run/mpl-cache" \
  /absolute/path/to/pymacruise/.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 inputs/t1.nii.gz --flair inputs/flair.nii.gz \
  --brain-mask inputs/mask.nii.gz --maseg inputs/maseg.nii.gz \
  --atlas-descriptor /absolute/path/to/pymacruise/tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir mini-run \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-convergence-threshold 0 --n-iter-max 2 --loop-thinning-limit 1000 \
  --finite-zero-variance
```

Append `--resume` with the same input paths and options to verify reuse. The
wheel was built from local commit `429eeb8`; the Python implementation is
unchanged from the full-suite-tested implementation commit
`c8a33033e823d9a2fbce5170ee7ebd984008e801`. That full suite reported
**1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,187.29
seconds**, peak sampled RSS **13.539 GB** under the 24 GB watchdog; see the
[full-suite report](full-pytest-2026-10-05-toads-ordering-carry.md).

## Configuration attempt and evidence boundary

An initial `workflow-run` config attempt used the same very small phantom with
the bundled full-size atlas, zero initialization-alignment iterations, and
reached a zero robust intensity span after atlas cropping. This was a
synthetic fixture/atlas-coverage mismatch, so the successful installed-wheel
run uses the compact committed atlas fixture and the explicit `plus-full-run`
atlas option. Config-driven `workflow-run` from the wheel, whole-workflow
Java/MATLAB equivalence, and real-MRI validation remain open.
