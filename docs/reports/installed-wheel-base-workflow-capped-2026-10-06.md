# Installed-wheel synthetic base workflow — capped TOADS run (2026-10-06)

## Result and scope

Built a wheel from the implementation tree at `a672c630a99386df9cb5a66f834368f3c9581aee`, installed it to an isolated target under `build/`, changed the working directory outside the checkout package, and ran `pymacruise.cli base-run` with `PYTHONPATH` containing only that installed target. The run resolved the VandyCRUISE atlas descriptor and payloads from the wheel and completed all six base stages: N4, first TOADS, Dura, second TOADS, CRUISE core, and MaSeg refinement.

All six manifest stages completed. The manifest contains 38 output records; all recorded hashes verified. Its 30 NIfTI outputs are finite and RAS. A second invocation from the same installed wheel verified and skipped all six stages. Peak sampled RSS was **6.426 GB**, below the 16 GB watchdog limit.

Each TOADS pass was deliberately capped at two iterations for this packaging integration check. The manifest correctly records `energy_threshold_reached=false`, `iteration_limit_reached=true`, and `converged=false` for both passes. This is installed-wheel stage, asset, output, and resume evidence; it does not establish TOADS convergence, Java/MATLAB equivalence, or real-data behavior. The separate source-checkout brain-resolution run remains the converged synthetic workflow evidence.

## Reproduction

From the repository root, build and install the wheel in the project-local ignored build directory:

```bash
mkdir -p build/wheel-base-workflow/dist build/wheel-base-workflow/site build/wheel-base-workflow/tmp
TMPDIR="$PWD/build/wheel-base-workflow/tmp" ./.venv/bin/python -m build --wheel --no-isolation --outdir build/wheel-base-workflow/dist .
./.venv/bin/python -m pip install --no-deps --target build/wheel-base-workflow/site build/wheel-base-workflow/dist/pymacruise-0.1.0-py3-none-any.whl
cd build/wheel-base-workflow
```

Run the installed package from this directory. The inputs are committed synthetic brain-resolution `.nii.gz` files; the atlas descriptor and referenced maps are read from the installed wheel:

```bash
PYTHONNOUSERSITE=1 PYTHONPATH="$PWD/site" MPLCONFIGDIR="$PWD/mpl-cache" \
  /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python \
  /Users/landmaba/dev/mipav/pymacruise/scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --atlas-descriptor "$PWD/site/pymacruise/data/atlases/cruise-atlas-10obj-toads2008.txt" \
  --maseg /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --output-dir run-capped \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 2 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Append the same options plus `--resume` to verify and skip the six completed stages. The recorded resume printed `Resume: verified <stage>; skipping` for each stage and exited successfully.

## Validation evidence

- Installed import path: `.../build/wheel-base-workflow/site/pymacruise/__init__.py`.
- Manifest: `status=complete`; six stages complete; both TOADS passes report the two-iteration cap and `converged=false`.
- Integrity: all 38 output records' SHA-256 hashes matched the manifest; 30 NIfTI outputs loaded with finite values and `('R','A','S')` axis codes.
- Resume: all six stages were hash-verified and skipped by the installed wheel.
- Memory: peak sampled RSS **6.426 GB** under a 16 GB watchdog; watchdog did not terminate the run.
- The implementation source did not change between full-suite commit `9f44477` and this wheel build. The latest full-suite evidence remains **1,207 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,217.29s**, peak sampled RSS **14.615 GB**, on `9f44477`; see [the full-suite report](full-pytest-2026-10-06-approx-thinning-validation-boundary.md).
