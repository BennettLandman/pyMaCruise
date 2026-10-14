# Human MaSeg surface-refinement compatibility comparison — 2026-10-12

> **Correction, 2026-10-14:** Fresh execution of the actual MATLAB routine shows that positive-integer `tabulate` includes zero-count rows, so `maxind` equals the modal label. The historical metrics below remain observations of the then-current Python implementations, but the claim that the low labels reproduce a MATLAB bug is incorrect. The Python default's distinct-label-row assumption was a port defect, now corrected. This explains the low-label discrepancy without requiring a different deployed MATLAB source, while remaining upstream differences still require comparison. See [fresh MATLAB evidence](matlab-refinement-direct-source-comparison-2026-10-14.md).

## Scope

The source-profile BASE run used the selected human inputs, identical level
sets, and `source_faithful=True`, which reproduces the historical MATLAB
`maxind` row-index assignment. I then reran only segmentation refinement on
the same MaSeg and level sets with `source_faithful=False`, assigning the
actual modal label instead. The rest of the run and all geometry are held
fixed. Both comparisons are by voxel index after exact input T1/MaSeg array
identity checks; they are not whole-workflow equivalence claims.

## Results

| Refinement behavior | Mean label Dice | Median | Minimum | Labels below 0.97 | Voxels differing from source | Python-only labels |
|---|---:|---:|---:|---:|---:|---|
| Historical `maxind` row-index behavior (`source_faithful=True`) | 0.91878 | 0.95992 | 0.00000 | 86/139 | 77,259 | 1, 2, 3, 5, 6 |
| Modal-label correction (`source_faithful=False`) | 0.96516 | 0.97551 | 0.38298 | 56/134 | 47,306 | none |

The corrected option reduces differing voxels by **29,953** and improves mean
label Dice by **0.04638**. It removes the five spurious low-valued labels, but
the refinement still fails the user's label-wise 0.97 threshold. Labels 208,
209, and 156 remain among the largest mismatches. This fixes one identified
source of disagreement without establishing full refinement equivalence.

## Source-provenance finding

All seven readable copies of `MAsegRefineBySurf.m` located in the supplied
MATLAB source tree and its source-import export have the same SHA-256
(`4c901b974dd8b7851758d87cb4bb892b19bf0427409e6dea8b4d0eda668e63e0`). Each
contains the historical assignment `labels(...)=maxind`; the intended lookup
`foundtable(maxind,1)` is commented out in both `fillGM` and `fillWM`. The
human input's MATLAB Compiler cache contains a binary `MCR` payload at the
cached `.m` path, not readable source, so it cannot confirm whether the
executed code differed. The downloaded V3 `SegRefine/SEG_refine.nii.gz`
contains none of the five low-valued labels produced by the Python
source-faithful run, while the supplied output tree has no MATLAB/Singularity
run log or manifest to identify the exact refinement code or any subsequent
label cleanup. The local readable V3 sources and observed V3 output therefore
still disagree on this behavior. Keep the source-faithful option and default
unchanged until the executed V3 source or a user-approved compatibility
choice resolves that provenance gap.

## Reproduction

The source-faithful result is part of
[`in-vivo-base-orientation-only-source-profile-2026-10-12.md`](in-vivo-base-orientation-only-source-profile-2026-10-12.md).
The corrected result was generated from that run's materialized input files and
unchanged inner/outer level sets:

```sh
RUN=build/invivo-base-orientation-only-source-profile-2026-10-12
.venv/bin/python -m pymacruise.cli segmentation-refine \
  --segmentation "$RUN/resources/orientation-only-inputs/MaSeg.nii.gz" \
  --inner-levelset "$RUN/stages/cruise-core/attempt-01/inner.nii.gz" \
  --outer-levelset "$RUN/stages/cruise-core/attempt-01/outer.nii.gz" \
  --output "$RUN/refinement-bugfix/SEG_refine.nii.gz" \
  --summary "$RUN/refinement-bugfix/summary.json" \
  --no-source-faithful
```

The held-fixed source comparison is in the ignored `build/` tree. The
reproducible comparison command uses
`scripts/compare_human_source_output.py --voxel-index-only` with the Python
corrected refinement as `--python-output` and the same core/preparation
directories as the source-profile run.
