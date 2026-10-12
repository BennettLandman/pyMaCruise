# Optional orthogonal RAS regridding: focused validation

Implementation commit: `689465d` (`Add default orthogonal RAS regridding to workflows`).

## Behavior

Both restartable full-workflow runners accept `regrid_inputs` (default `true`),
and their CLIs expose `--regrid-inputs` and `--no-regrid-inputs`. Workflow JSON
accepts the same setting under the base run options or the plus `run_options.plus`
group. Axis-aligned RAS inputs pass through without resampling. For oblique RAS
inputs, T1 is resampled linearly to an orthogonal grid at its original spacing;
the matching user mask and MaSeg are resampled to the same grid with nearest
neighbor. The plus workflow leaves FLAIR on its own grid for the existing
registration step. Derived images are saved as `.nii.gz`, and their paths and
SHA-256 hashes are recorded in the top-level run manifest.

With regridding disabled, an oblique but orthogonal RAS T1/mask/MaSeg grid is
preserved without resampling. The selected TOADS registration adapter uses the
source MIPAV axis signs and voxel spacings; it does not pass affine rotation
terms into atlas registration. The full RAS affine remains on saved images for
physical-space operations and surfaces. Sheared grids remain unsupported.
Inputs with non-RAS axis codes are still reoriented to RAS without
interpolation and emit the standard warning.

## Evidence

Focused command:

```text
./.venv/bin/pytest -q tests/test_image_io.py tests/test_base_workflow_restartable.py tests/test_plus_full_workflow.py
17 passed in 118.48s
```

The focused set includes an oblique, restartable synthetic base workflow, shared
RAS-grid and nearest-neighbor checks, plus-workflow resource/hash checks, CLI
and workflow-config defaults, and early opt-out rejection. The two plus
regridding selectors and the base CLI default/opt-out selector were rerun on
the current checkout: **3 passed in 2.21s**. `pymacruise plus-full-run --help`
lists both regridding switches. `git diff --check` passed.

This is synthetic integration evidence. It does not establish source
equivalence or real-MRI output equivalence. The earlier full suite is recorded
separately in [the full-suite report](full-pytest-2026-10-11-orthogonal-regridding.md)
for tested tree `cabb52aa`; that run predates the current implementation tree
and is not a fresh full-suite result for the current checkout.

## Oblique no-resampling opt-out — 2026-10-12

The previous opt-out rejected oblique RAS grids. It now preserves an oblique,
orthogonal RAS input grid in both full workflows, records
`regrid_inputs=false`, and creates no regridded input resources. The TOADS
geometry adapter accepts rigid obliquity and still rejects shears. CLI help
describes the opt-out as preserving the original grid without interpolation.

Focused validation result: **19 passed, 0 failed, 0 skipped in 0.94s**, peak
sampled RSS **0.125 GB** under an 8 GB watchdog. The selected tests cover grid
classification, TOADS signed-spacing behavior, no-regrid manifest/resource
behavior for both workflows, default regridding, and shear rejection. This
does not yet compare a no-regrid oblique workflow with a complete historical
Java/MATLAB run.
