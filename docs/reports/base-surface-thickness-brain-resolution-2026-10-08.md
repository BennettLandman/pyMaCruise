# Base scientific surfaces and thickness on brain-resolution synthetic inputs — 2026-10-08

## Scope

The base V3 runner produces scientific inner/central/outer surfaces and thickness products. This report records the new `base-surface-thickness` stage adapter run on the full-grid level sets from the existing synthetic brain-resolution base run. The focused restart regression separately exercises the newly integrated seven-stage workflow on a smaller synthetic fixture. The previous brain-resolution run manifest predates this stage, so this direct full-grid command does not claim that the seven-stage brain-resolution manifest has been run end to end.

## Reproduction

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-surface-thickness \
  --inner-levelset build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-04/inner.nii.gz \
  --central-levelset build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-04/central.nii.gz \
  --outer-levelset build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-04/outer.nii.gz \
  --output-dir build/ace-indexed-heap-2026-10-08/base-surface-thickness \
  --provenance build/ace-indexed-heap-2026-10-08/base-surface-thickness/provenance.json
```

The command exited 0 and peaked at **0.577 GB RSS** under the 4 GB watchdog. Input level sets are `(172,216,164)` and their SHA-256 hashes are `31233e72558c488a7e536373571c12293f9ca64d464e84f324c0f4cff056d042` (inner), `4fdaf13b1d1ffd80df8e7c0e2a6c3b2290cc761da16148987166ecfa6aa0445c` (central), and `1f88ba7b4483fa6e59af759055db54577046cefcb5516045e5b7d469a3ace6c2` (outer).

## Outputs and checks

| Product | Result |
|---|---|
| Inner surface | RAS-mm VTK; 116,178 vertices, 232,368 faces |
| Central surface | RAS-mm VTK; 121,100 vertices, 242,224 faces |
| Outer surface | RAS-mm VTK; 123,872 vertices, 247,828 faces |
| Inner surface with point data | RAS-mm VTK; stores inner-to-outer distance, return distance, and their mean at each of 116,178 vertices |
| Level-set thickness | finite float32 RAS `.nii.gz` on `(172,216,164)` grid; range 0–5.2573075 mm; 142,858 crossing voxels |
| Summary | Level-set crossing mean 2.00259464 mm (sample SD 0.35367777 mm); closest-point surface mean 1.64906598 mm (sample SD 0.49805758 mm) |

The measurements are synthetic algorithm outputs; no thickness-specific owner acceptance threshold exists. Surface extraction uses scikit-image Lewiner marching cubes and does not reproduce JIST `IsoSurfaceOnGrid` connectivity 18. The held-fixed JIST thickness wrapper comparisons do not eliminate this mesh-connectivity risk. No renderer or viewer product is created.

## Focused validation

```sh
./.venv/bin/python -m pytest -q tests/test_base_workflow_restartable.py
```

Result: **3 passed, 0 failed, 0 skipped in 34.48s**. The test runs the seven-stage workflow on a smaller synthetic fixture, checks RAS/finite volume outputs and RAS-mm meshes with thickness point data, verifies all stage checkpoints, confirms same-settings resume, and checks downstream invalidation after core-output corruption. No full pytest suite was run in this cycle.

The previous brain-resolution six-stage run and its resume evidence are recorded in [the 2026-10-08 restart audit](base-workflow-restart-after-thickness-change-2026-10-08.md). Full-size seven-stage manifest validation, JIST mesh extraction equivalence, whole-workflow source parity, plus parity, and real-MRI validation remain open.
