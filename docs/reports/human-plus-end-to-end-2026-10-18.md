# Human-data PLUS end-to-end run — 2026-10-18

## Result

The Python PLUS pipeline completed end to end on the supplied human case
100260. The top-level `plus-full-workflow.json` and all eleven configured
stages report `complete`. The run produced the final refined MaSeg,
three full-grid level sets, inner/central/outer VTK surfaces, and a thickness
volume. A second invocation with the same inputs and `--resume` returned the
same completed outputs without recomputing the expensive stages.

This is human-execution evidence for the current Python workflow. It is not an
independent execution of the original MaCRUISE source and does not establish
full source equivalence or clinical validity.

## Inputs and selected atlas

Inputs were the supplied T1, FLAIR, brain mask, and MaSeg under
`source-reference-run/inputs/100260/`. The AtlasEM descriptor was the adapted
`ibsr11objV3` atlas at
`build/ibsr11objV3-python-ras/LesionToads-11-obj-atlasV3-atlasem-python.txt`;
the workflow also used the 12-object LesionTOADS descriptor and its topology
image. Input hashes are recorded in
`build/plus-human-100260-atlasem-2026-10-18/run-numba-v1/plus-full-workflow.json`.

The final T1/RAS grid is 176 × 240 × 256 with 1 mm spacing. Final MaSeg,
level sets, and thickness volume all have that shape and affine. The surface
workflow's crop was 140 × 178 × 134.

## Execution evidence

- AtlasEM and Dura initialization completed with covariance floor `1e-6`;
  the floor regularizes singular class covariances encountered on this case.
- The registered lesion loop stopped by convergence on iteration 26, with
  energy ratio `5.6365363e-5` below the configured `1e-4` threshold; it did
  not stop at its iteration limit.
- Selected topology, ACE, GVF, and nested TGDM completed, followed by
  uncropping, surfaces, thickness, and final MaSeg refinement.
- All emitted NIfTI arrays checked here contain finite values. The final
  segmentation is integer-valued and contains the observed MaSeg labels from
  0 through 209. All three surfaces are nonempty meshes (about 209k–215k
  vertices each).
- Thickness volume: 0–10.244 mm, mean 2.459 mm over crossing voxels;
  surface thickness summary: mean 1.236 mm over 210,174 vertices. These are
  descriptive values, not source-comparison acceptance results.
- Focused regression tests for AtlasEM statistics/iterations and lesion-shape
  registration: **25 passed**.

## Reproduction

From the `pymacruise` repository root:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mpl .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 source-reference-run/inputs/100260/T1.nii.gz \
  --flair source-reference-run/inputs/100260/FLAIR.nii.gz \
  --brain-mask source-reference-run/inputs/100260/brain_mask.nii.gz \
  --maseg source-reference-run/inputs/100260/MaSeg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --atlasem-descriptor build/ibsr11objV3-python-ras/LesionToads-11-obj-atlasV3-atlasem-python.txt \
  --no-regrid-inputs \
  --output-dir build/plus-human-100260-atlasem-2026-10-18/run-numba-v1 \
  --resume
```

The primary run manifest is
`build/plus-human-100260-atlasem-2026-10-18/run-numba-v1/plus-full-workflow.json`.
Surface products and per-stage hashes are under
`build/plus-human-100260-atlasem-2026-10-18/run-numba-v1/code-revisions/74bd4c6b9f79dc3d/plus-surface/stages/`.

## Remaining limits

The FLAIR input has oblique RAS axes. Because this run used `--no-regrid-inputs`,
the workflow reoriented voxel axes without interpolation and discarded those
oblique direction cosines; FLAIR voxel-center/world correspondence is not
preserved. The run emits a warning for this. A human run intended for
anatomical use should be repeated with an explicit, physically correct FLAIR
registration/regridding path.

This run has not been compared end to end against original MaCRUISE on the
same human inputs. The Python registration path is not verified as numerically
equivalent to the source Aladin optimizer, AtlasEM's covariance floor changes
singular-matrix behavior, and the selected atlas adapter has not been checked
against a complete original-source AtlasEM run. Surface and thickness values
above therefore establish successful execution and valid output structure,
not equivalence or accuracy.
