# Installed-wheel synthetic MaCRUISEplus workflow (2026-10-08)

## Result and scope

Built `pymacruise-0.1.0-py3-none-any.whl` from commit `9e7484ff9e19819179cf5e2223f7abe1ff97cd84`, installed it with `pip --no-deps` into an isolated target, and ran the workflow from a build directory outside the source checkout. `PYTHONNOUSERSITE=1` and `PYTHONPATH` pointed only to the installed package. The mini LesionTOADS atlas descriptor and maps were supplied from the committed synthetic fixture; package imports and atlas resources resolved successfully.

The raw-input plus workflow finished all nine stages and reported `status=complete`, `terminal=true`. LesionTOADS converged at iteration **6**, with energy ratio **`7.722568398094154e-07`** below the configured `1e-4` threshold; it did not stop at its iteration limit. Final products include native-T1 tissue transfers, three RAS-mm surfaces, level sets, and thickness volume/metrics.

An independent recursive manifest audit verified **86 unique manifest-referenced file hashes**. All **62 NIfTIs** were finite and RAS. All three VTK meshes had finite vertices, triangular faces, and valid point indices. The complete run peaked at **0.729 GB RSS** under a 4 GB watchdog. Same-settings `--resume` reported the raw-input workflow verified/resumed, exited 0, and peaked at **0.234 GB RSS**.

This is isolated-wheel packaging, bundled-resource, synthetic workflow, output-integrity, natural-convergence, and restart evidence. It does not establish composite Java/MATLAB equivalence, cross-grid surface-position parity, or real-MRI behavior.

Wheel SHA-256: `0ecc7029083c70d095d97ebe16b9308c53cd67cbfc7f89cc951800bfa97c8c7a`.

## Reproduction

From the repository root, build and install the wheel:

```bash
PROJECT=/Users/landmaba/dev/mipav/pymacruise
mkdir -p "$PROJECT/build/wheel-plus-workflow-current/dist" \
  "$PROJECT/build/wheel-plus-workflow-current/site" \
  "$PROJECT/build/wheel-plus-workflow-current/tmp"
TMPDIR="$PROJECT/build/wheel-plus-workflow-current/tmp" \
  "$PROJECT/.venv/bin/python" -m build --wheel --no-isolation \
  --outdir "$PROJECT/build/wheel-plus-workflow-current/dist" "$PROJECT"
"$PROJECT/.venv/bin/python" -m pip install --no-deps --target \
  "$PROJECT/build/wheel-plus-workflow-current/site" \
  "$PROJECT/build/wheel-plus-workflow-current/dist/pymacruise-0.1.0-py3-none-any.whl"
cd "$PROJECT/build/wheel-plus-workflow-current"
```

Run outside the checkout. The T1, FLAIR, brain mask, and MaSeg are committed in `tests/data/synthetic/plus_scaled_grid/`; the mini atlas descriptor is in `tests/data/synthetic/plus_lesion_atlas_init/`:

```bash
PROJECT=/Users/landmaba/dev/mipav/pymacruise
FIXTURE="$PROJECT/tests/data/synthetic/plus_scaled_grid"
ATLAS="$PROJECT/tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt"
env PYTHONNOUSERSITE=1 PYTHONPATH="$PWD/site" MPLCONFIGDIR="$PWD/mpl-cache" \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$PROJECT/.venv/bin/python" "$PROJECT/scripts/run_with_rss_watchdog.py" \
  --limit-gb 4 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  "$PROJECT/.venv/bin/python" -m pymacruise.cli plus-full-run \
  --t1 "$FIXTURE/t1.nii.gz" --flair "$FIXTURE/flair.nii.gz" \
  --brain-mask "$FIXTURE/brain-mask.nii.gz" --maseg "$FIXTURE/maseg.nii.gz" \
  --atlas-descriptor "$ATLAS" --output-dir run-wheel-validated \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

Repeat the same command with `--resume` and a 4 GB watchdog to verify the top-level workflow restart.

## Interpretation

- The run used only the installed wheel on `PYTHONPATH`; imaging dependencies came from the project virtual environment.
- The plus-specific synthetic workflow naturally converged and wrote its full nine-stage product set.
- The output audit covers manifests and intermediate/final products under this run root; duplicate references to the same file were counted once.
- Selected component comparisons and known deviations remain documented separately. Synthetic completion must not be read as full Java/MATLAB equivalence.
