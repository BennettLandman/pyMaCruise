# Brain-resolution synthetic base workflow with full registration — 2026-10-08

## Reproduction

This run exercises the complete seven-stage base workflow on the committed brain-resolution synthetic T1, matching user-provided brain mask, and MaSeg fixture. It uses all four TOADS registration-pyramid levels and permits up to 50 iterations in each TOADS pass and registration phase.

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/base-workflow-2026-10-08-natural-registration \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 50 --registration-max-iterations 50 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

The run completed in **3:08:33**. All seven stages passed, both TOADS passes met their energy stopping criteria before the 50-iteration cap, and the manifest reports `status=complete`, `converged=true`. The first pass converged at iteration 11 (`8.55285762e-05`); the second converged at iteration 48 (`5.68501309e-05`). The 16 GB RSS watchdog sampled a peak of **8.514 GB** and did not terminate the run.

## Outputs and restart

An independent audit checked **46/46 manifest output hashes and byte counts**, including **31 finite RAS NIfTIs**. All seven stage records are complete; the manifest lists 24 final products, including CRUISE core level sets and ACE/GVF products, the refined segmentation, and surface/thickness outputs. The inner, central, and outer meshes contain 115,516, 120,810, and 123,944 vertices. Level-set thickness has 142,549 crossing voxels and mean 1.99991 mm; symmetric closest-point surface thickness has mean 1.65019 mm.

The same command and options with `--resume` exited 0 and verified/skipped all seven stages; peak sampled RSS was **0.121 GB**. The per-file hashes, NIfTI shapes/orientations, input hashes, options, and stage stop metadata are in [machine-readable evidence](base-workflow-natural-registration-2026-10-08.json). The generated run directory is ignored build output and is reproducible from the command above.

## Selected Java source comparisons on this run

The saved full-registration outputs were compared against the supplied Java sources on held-fixed inputs:

- **TGDM core:** `TGDMCoreReferenceHarness` ran the selected `MaGenericTGDM`/`GenericTGDM` calls on the exact six RAS crop inputs emitted by this run, using the selected one-iteration wrapper profile. On the native `(172,216,164)` grid at `0.9375 mm` spacing, inner, central, and outer Dice were `1.0`, `1.0`, and `0.9999894253`. Symmetric mean nearest-surface distances were `2.67e-06`, `1.32e-06`, and `0.0001399 mm`; Hausdorff distances were `0.00196`, `0.00329`, and `0.08876 mm`. Each output passes the agreed Dice ≥0.97, symmetric mean distance ≤0.01 mm, and Hausdorff ≤1 mm gates. The comparison covers the delivered masks, not every intermediate evolution field.
- **JIST surfaces:** `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0.0f, false)` consumed the exact saved inner, central, and outer fields. Each Java/Python mesh had identical vertex and face counts, bijective vertex correspondence, and 100% ordered triangle matches. The maximum paired-vertex displacement was `7.1412e-06 mm` for each surface; both meshes for each surface had zero boundary or nonmanifold edges.

Reproduce the surface comparison with:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py \
  --core-dir build/base-workflow-2026-10-08-natural-registration/stages/cruise-core/attempt-01 \
  --work-dir build/base-workflow-natural-registration-jist18-2026-10-08
```

The TGDM comparison was reproduced by preparing the six crop inputs with `scripts/prepare_tgdm_core_java_inputs.py`, compiling `tests/java_reference/TGDMCoreReferenceHarness.java` with the stubs and selected upstream source under `tests/java_reference`, running the harness with a 4 GB RSS watchdog (`-Xmx3g`), then running `scripts/convert_tgdm_core_java_reference.py` and `scripts/compare_tgdm_core_saved_java_reference.py` against this run's `cruise-core` directory. The Java harness exited 0 and peaked at **1.754 GB** sampled RSS. Full per-output metrics, exact field hashes, Java source/harness hashes, and reproducibility commands are in the machine-readable evidence.

## Interpretation and limits

This establishes a restartable, naturally converged, full-registration synthetic base workflow using the bundled 11-class TOADS atlas, with selected TGDM outputs and JIST surface extraction compared on this exact run. It does not establish whole-workflow MATLAB/JIST equivalence. The base wrapper deliberately used one Dura pass, one ACE diffusion iteration, 200 GVF iterations, and one TGDM iteration for this bounded end-to-end run; these are not claims that every source default was exercised. Whole-composite CRUISE/refinement parity and real-MRI validation remain open. The supplied human T1 still needs its matching brain mask.

The most recent full pytest suite on the same implementation code (tested tree `5a0adcf`, before documentation-only commits) passed **1,502 tests, 0 failed, 0 skipped, and 11 expected RAS warnings in 1,873.83 s**, with **12.442 GB** peak sampled RSS. No code changed after that suite; this workflow run adds no newly mapped selected VandyCRUISE TOADS source lines.
