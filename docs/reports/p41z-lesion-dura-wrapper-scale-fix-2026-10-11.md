# P41z LesionTOADS Dura wrapper scale correction — 2026-10-11

## Source finding and implementation

The selected Java `MedicAlgorithmLesionToads` plugin saves the CSF, GM,
white-matter fill, and white-matter mask products after calling
`ImageDataMath.scaleFloatValue(..., 255.0f)` on each channel (lines 509–512 of
the supplied source). The low-level `generateDuraRemovalOutputs` method
returns normalized channel values and its Java harness comparison remains
valid for that method. The Python checkpoint adapter had omitted the plugin
wrapper's subsequent scale, so its products were not the images consumed by
the source MaCRUISEplus pipeline.

`run_plus_lesion_dura_output_checkpoint` now leaves the low-level generator
unchanged and applies a float32 ×255 scale at the wrapper boundary. The
checkpoint signature and source-order manifest record this operation so
previous unscaled products cannot be silently resumed as valid results.
The end-to-end synthetic workflow regression checks all four saved products
against the raw generator outputs multiplied by 255.

## Validation

Focused checks passed: **3 passed in 39.93 seconds**, covering the complete
restartable synthetic plus workflow, the regridding defaults/opt-out, and the
base regridding CLI defaults/opt-out.

The correction was also applied to a held-fixed P41z iteration-11 checkpoint
and passed through the downstream plus surface workflow. The source iteration
remained nonterminal; only Dura wrapper products and downstream stages were
recomputed. MaSeg label support changed as follows:

| Check | Before wrapper scale | After wrapper scale | Acceptance |
| --- | ---: | ---: | ---: |
| WM label 44 inside inner surface | 59.51% | 99.406% | ≥97% |
| GM label 100 between inner and outer surfaces | 0.112% | 92.835% | ≥97% |
| CSF label 5 outside outer surface | 99.963% | 99.258% | ≥97% |

Level-set subtraction thickness mean improved from **0.909 mm** to **4.832
mm**, within the selected 5±1 mm range. The new meshes are closed, but each
has two connected components and Euler characteristic −12. Surface nearest
distance median is **6.434 mm** (5th–95th percentile 0.237–12.280 mm), outside
the selected median and percentile bounds. Thus the wrapper correction
explains most of the prior WM-support loss and restores plausible mean
thickness, but this checkpoint still fails ribbon acceptance because GM
support, mesh topology, and surface thickness bounds remain out of range.

Machine-readable metrics and outputs are under the ignored `build/`
directory at
`build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/`.
This is a held-fixed synthetic checkpoint integration, not an independently
generated Java workflow, whole-loop parity, or real-human-MRI validation.

## Regridding option

Both `base-run` and `plus-full-run` already support default-on, recommended
orthogonal RAS regridding and the explicit `--no-regrid-inputs` option. RAS
inputs already on an axis-aligned grid are not interpolated. Disabling
regridding currently requires an axis-aligned RAS T1 because the Python TOADS
registration adapter does not yet accept oblique/sheared grids. Therefore,
the opt-out reproduces an earlier axis-aligned execution environment; it does
not yet enable the original oblique human input to run without regridding.
