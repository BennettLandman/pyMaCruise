# Corrected-code synthetic MaCRUISEplus rerun

## Outcome

The folded 64×80×64 synthetic plus workflow completed all nine stages on
implementation commit `b785bff` (`Match retained TOADS registration center`).
The registered LesionTOADS loop naturally stopped at iteration 6 with energy
ratio `7.722568398094154e-07`; the iteration limit was not reached. The top
manifest records `status=complete` and `terminal=true`.

An independent audit verified all 92 unique file paths with SHA-256 records
referenced by the nested manifests (including inputs and outputs). All 67
referenced NIfTIs were finite and RAS-oriented. The three algorithmic VTK
surface products declare RAS coordinates in millimetres. A same-settings
`--resume` verified and resumed the completed workflow successfully under a
4 GB watchdog, peaking at 0.237 GB. The initial run peaked at 0.800 GB under a
16 GB limit and took approximately 51 minutes from output-directory creation
to final manifest completion; most of the elapsed time was the registered
LesionTOADS loop, which currently provides no per-iteration checkpoint while
it runs.

This confirms synthetic implementation and restart behavior after the shared
TOADS retained-center correction. It does not establish whole-workflow
Java/MATLAB equivalence or real-MRI validation. Per-iteration progress and
runtime optimization for the registered LesionTOADS loop remain open.

## Reproduction

```sh
env PYTHONNOUSERSITE=1 PYTHONPATH=src \
  MPLCONFIGDIR=build/plus-retained-center-2026-10-09/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-retained-center-2026-10-09 \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

The resume command repeats the same arguments and appends `--resume`, using a
4 GB watchdog. It reports `Verified and resumed raw-input plus LesionTOADS
workflow` and exits successfully.

The full repository suite on the same implementation commit passed **1,503
tests, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,764.57
s**, with **13.787 GB** peak sampled RSS; see [full-suite evidence](full-pytest-2026-10-09-retained-center.md).

## Held-fixed selected-source comparisons on this run

The run's exact saved crop and level sets were replayed through selected Java
or JIST source harnesses. These checks isolate the named stages; they do not
compare the full upstream pipeline end to end.

| Behavior | Result on the current run |
|---|---|
| Plus tissue fusion (`CombineMAandTOADS`) | GM, WM, and mask are exact across all 331,056 crop voxels. |
| Plus topology correction and MaSeg refinement | Topology-corrected image, MaWM erosion, and MaSeg-refined output are exact across all 331,056 voxels. |
| Plus raw atlas boundary (`MinFilter3D`) | The saved crop is empty on both sides (0 changed voxels; vacuous). The separate folded source fixture is nonempty (17,190 foreground voxels) and exact. |
| Plus ACE mode 3 | Enhanced GM, skeleton, and thinned skeleton are exact across all 331,056 voxels. |
| Plus GVF | Maximum absolute vector-component error `2.1234e-7`; mean `2.3023e-8`, within existing Java-reference limits. Java peak RSS was 0.247 GB. |
| Plus TGDM | Central/inner/outer Dice is 1.0. Symmetric mean surface distances are `1.6512e-6`, `2.1938e-7`, and `1.6877e-6` mm; Hausdorff distances are `0.0003467`, `0.0000866`, and `0.0003281` mm. Java peak RSS was 0.358 GB. |
| JIST `IsoSurfaceOnGrid` versus Python surfaces | Each of the three surfaces has 3,724 vertices and 7,444 faces; every ordered face matches, topology matches (one component, Euler characteristic 2), and maximum paired-vertex distance is `8.72e-7` mm. |

The checked-in reproductions are `compare_plus_tissue_update_saved_java.sh`,
`compare_plus_topology_saved_java.sh`, `compare_plus_boundary_saved_java.sh`,
`compare_plus_ace_gvf_saved_java.sh`, `compare_plus_tgdm_saved_java.sh`, and
`compare_workflow_jist_surfaces.py`; each consumed the run under
`build/plus-retained-center-2026-10-09/`. The Java evidence holds upstream
inputs fixed at the saved Python stage products. It therefore strengthens
stage-level source validation, but does not establish LesionTOADS or complete
plus-workflow Java/MATLAB equivalence.
