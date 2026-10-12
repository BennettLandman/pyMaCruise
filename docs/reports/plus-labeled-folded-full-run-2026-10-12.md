# Labeled folded synthetic PLUS full-run attempt

## Result

The current `plus-full-run` path completed input handling, orthogonal RAS
regridding, preprocessing, registration, registered recropping, first
membership initialization, and first thinning on the committed 65³ folded
synthetic fixture. It then stopped during the post-thinning GM-enabled
membership update because the source-faithful objective became nonfinite.
The manifest correctly records the run as failed. This is not a successful
end-to-end PLUS run and does not establish ribbon quality or source equivalence.

The fixture has labels 0/5/44/100 and deterministic tissue-specific T1/FLAIR
values, but it contains no lesion example. The bundled atlas descriptor
enables lesion handling based on its topology declarations. The V3 MATLAB
wrapper also explicitly selects `With Lesion`, and VandyCRUISE sets
`hasLesion=true` from the lesion atlas topology; absence of a lesion seed does
not select the no-lesion atlas. The first-iteration summary records two zero
lesion-standard-deviation channels and no zero standard-deviation tissue
channels. `computeVariances` sets an empty lesion variance to zero, and the
existing Java selected-method zero-variance harness confirms the resulting
similarity can be NaN. Python's source-faithful numerical behavior therefore
exposes a nonfinite objective, while this Python workflow currently stops at
that point instead of continuing the Java loop. That stop-vs-continue behavior
is an explicit integration gap; NaNs must not be silently replaced. A separate
lesion-bearing fixture can exercise finite loop behavior, while source-faithful
handling needs direct continuation/output comparison. The finite-zero-variance
switch remains an explicit comparison option and was not used here.

Source anchors: `run_LesionToads_jist.m:37` selects `With Lesion`;
`LesionToadSegmentation.java:230-252` derives `hasLesion` from atlas topology;
`computeVariances` sets an empty lesion standard deviation to zero at lines
2738-2741. The existing selected-method zero-standard-deviation fixture at
`tests/data/synthetic/lesion_toads_similarity/java_reference.json` records the
Java NaN case.

## Resources and artifacts

Peak sampled RSS was 3.378 GB under the 16 GB watchdog, so the failure was
algorithmic rather than a memory limit. Completed stages and the failure
manifest are retained under the ignored `build/full-plus-labeled-folded-current/`
directory. In particular, see:

- `plus-full-workflow.json`
- `lesion_workflow/first-iteration/first_iteration_summary.json`
- `lesion_workflow/registered-branch/lesion-registered-branch.json`

The input fixture remains at
`tests/data/synthetic/plus_labeled_folded_cortex/`. This run is synthetic
integration evidence only. It is distinct from the smaller CLI integration
test that completed but failed its independent ribbon-quality audit, and from
all Java source-comparison and real-human evidence.

## Reproduction

From the repository root:

```bash
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_labeled_folded_cortex/t1.nii.gz \
  --flair tests/data/synthetic/plus_labeled_folded_cortex/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_labeled_folded_cortex/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_labeled_folded_cortex/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/full-plus-labeled-folded-current
```

## Next step

Implement a source-anchored Java comparison for the zero-variance continuation
and its exported maps, then create a lesion-bearing synthetic case to exercise
finite full-loop behavior. Keep source-faithful NaN behavior enabled and make
the Python workflow's stop-versus-continue policy explicit until the Java
comparison establishes the expected outputs.

## Stabilized one-pass follow-up

To determine whether the complete workflow wiring can finish beyond this
source-faithful stop, the same inputs were run with the explicit
`--finite-zero-variance` option and `--max-passes 1`. This is a deliberately
stabilized, non-source-faithful comparison and is not used to claim Java
equivalence.

```bash
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_labeled_folded_cortex/t1.nii.gz \
  --flair tests/data/synthetic/plus_labeled_folded_cortex/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_labeled_folded_cortex/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_labeled_folded_cortex/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/full-plus-labeled-folded-finite-onepass \
  --max-passes 1 --finite-zero-variance
```

All nine public workflow stages completed and wrote the expected LesionTOADS,
Dura, level-set, mesh, and thickness products. The run is complete but
nonterminal at iteration 2 (energy ratio `0.033323`, threshold `0.0001`). Peak
sampled RSS was 5.046 GB under the 16 GB watchdog. An independent audit found
all 95 manifest hash records valid and all 64 NIfTIs finite on RAS grids.
Ribbon quality **failed**: GM between inner and outer support is 11.28% (target
≥97%); mean level-set thickness is 3.370 mm (target 4–6 mm); meshes are
watertight but have 13 connected components and Euler characteristic 10. The
surface closest-point median/p05/p95 is 0.173/0.140/2.253 mm. The audit report
is retained at `build/full-plus-labeled-folded-finite-onepass/ribbon-quality.json`.

The main performance cost is also measured: the iteration-2 ordered membership
sweep took 1,014 seconds, while thinning took 42 seconds, statistics refresh
39 seconds, and shape registration/branching 323 seconds. This result
establishes a complete bounded stabilized workflow and a concrete runtime
target, but it fails anatomy criteria and is not source-faithful, converged,
Java-equivalent, or real-human validation.
