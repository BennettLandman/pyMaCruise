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

With regridding disabled, an oblique T1 is rejected before a run manifest is
created because the current TOADS registration implementation requires an
axis-aligned RAS grid. An already axis-aligned RAS input can be run with
`--no-regrid-inputs` without interpolation.

## Evidence

Focused command:

```text
./.venv/bin/pytest -q tests/test_image_io.py tests/test_base_workflow_restartable.py tests/test_plus_full_workflow.py
17 passed in 118.48s
```

The focused set includes an oblique, restartable synthetic base workflow, shared
RAS-grid and nearest-neighbor checks, plus-workflow resource/hash checks, CLI
and workflow-config defaults, and early opt-out rejection. The two newly added
plus regridding selectors were rerun after strengthening the manifest-hash
assertion: **2 passed in 1.08s**. `pymacruise plus-full-run --help` lists both
regridding switches. `git diff --check` passed.

This is synthetic integration evidence. It does not establish source
equivalence or real-MRI output equivalence. The post-change full pytest suite
has not yet run; roadmap item P46 remains Started until that run is recorded.
