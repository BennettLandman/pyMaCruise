# N4 no-mask rescaling profile correction — 2026-10-14

## Finding

The MaCRUISE V3 MATLAB wrapper calls `N4BiasFieldCorrection -i input -o output` without `-x` (mask) or `-r` (rescale). The ANTs command-line implementation initializes rescaling to true, but then sets `doRescale=false` whenever no mask was explicitly specified. The command's option documentation also says a mask is required for rescaling. Therefore the selected wrapper profile has no mask and rescaling **off**.

The ANTs source branch can be inspected in [`N4BiasFieldCorrection.cxx`](https://github.com/ANTsX/ANTs/blob/main/Examples/N4BiasFieldCorrection.cxx); this records the present upstream logic. The selected MATLAB wrapper invocation is documented in [the source-output comparison](in-vivo-source-output-comparison-2026-10-11.md) and source map. Exact version of the bundled historical executable remains unavailable.

## Code correction

`N4Config` and the standalone `n4` CLI now default `rescale_intensities=False`. An explicit `--rescale-intensities` remains available for profiles that intentionally request the post-step. The default no-mask call preserves the wrapper's source command behavior. This corrects the former inference that ANTs CLI's internal `doRescale=true` meant the no-mask invocation rescaled its output.

## Spline-profile check

The bare ANTs CLI's documented B-spline default is one mesh element per image
axis (`-b [1x1x1,3]`). ANTsPy documents `spline_param=None` as that same
one-element mesh in every dimension. The installed ANTsPy 0.6.3 adapter
confirms this directly: `None` becomes `[1] * image.dimension`, and a
three-element list is serialized as a per-axis mesh. PyMaCRUISE explicitly
passes `[1, 1, 1]`, so this is equivalent to the selected bare-CLI mesh; this
parameter does not explain the remaining historical N4 difference. The
historical ANTs build is unavailable, so equivalence to that exact binary
cannot be asserted.

Sources: [ANTs CLI N4 source and option defaults](https://github.com/ANTsX/ANTs/blob/main/Examples/N4BiasFieldCorrection.cxx),
[ANTsPy N4 API](https://antspy.readthedocs.io/en/latest/core.html), and the
installed `ants/ops/bias_correction.py` implementation (ANTsPy 0.6.3).

## Validation

Focused command:

```bash
MPLCONFIGDIR=build/test-mpl PYTHONNOUSERSITE=1 PYTHONPATH=src \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/pytest -q tests/test_n4.py tests/test_cli_n4.py \
  tests/test_plus_lesion_preprocessing.py
```

Result: **8 passed, 0 failed, 0 skipped in 7.76s**, one expected RAS reorientation warning. This is focused adapter/CLI evidence, not a full-suite result. The no-mask human N4 CLI completed on the supplied `(184,256,256)` T1 using **0.669 GB peak sampled RSS** under an 8 GB watchdog. Compared with the historical N4 volume after axis-only RAS handling, MAE is **6,301.90**, RMSE **14,898.90**, and maximum absolute difference **138,980.28** intensity units; historical and Python ranges are `[0, 210566.03]` and `[0, 71585.75]`. This does not match the historical intermediate and does not establish source equivalence; the original run configuration/build remains unknown. In particular, the previous rescale-on human comparison was an exploratory measurement, not justification for the bare no-mask wrapper's defaults.

The corrected N4 output is saved only under ignored `build/invivo-n4-source-profile-correction-2026-10-14/`. The prior complete human BASE run used rescaling on, so its workflow outputs remain valid evidence for that explicit profile but are no longer the current selected wrapper-profile integration. A fresh full BASE run with rescaling off is the next integration step.
