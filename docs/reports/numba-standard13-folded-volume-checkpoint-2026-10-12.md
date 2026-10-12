# Ordered standard-13 membership acceleration and PLUS checkpoint

## Result

The standard-13 TOADS membership sweep now has a Numba implementation and keeps the prior Python implementation as an executable oracle (`backend="python"`). Numba is the default. The kernel preserves source traversal and in-place neighbor updates, class and energy state ordering, float32 arithmetic, the six-face neighborhood, and the source-selected zero-variance behavior.

The focused oracle comparison is bit-for-bit on the 5³ Java-harness input and on a 15³ folded volume tiled from that fixture. All returned membership, relation, energy, denominator, objective-history arrays and the scalar objective match exactly. The focused regression group passed **11 tests in 6.49 s**. This is Python-backend equivalence evidence on synthetic inputs; it is not a new Java-equivalence claim.

## Integrated synthetic PLUS run

Command (with an RSS watchdog capped at 16 GB):

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 30 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_labeled_folded_cortex/t1.nii.gz \
  --flair tests/data/synthetic/plus_labeled_folded_cortex/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_labeled_folded_cortex/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_labeled_folded_cortex/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-labeled-folded-numba-onepass-2026-10-12 \
  --max-passes 1
```

All nine public workflow stages completed, with peak sampled RSS **5.327 GB**. The run is a deliberate one-pass checkpoint (`terminal=false`), not a converged result. In the iteration-2 progress manifest, `update_memberships` took **1.10 s**. An earlier run of the same labeled-folded 65³ case recorded about **1,014 s** for that stage. Treat the ratio as an indicative performance measurement, not a controlled benchmark: comparison of saved outputs found that upstream registered FLAIR already differs (maximum absolute difference 0.65204) even though both run manifests list the same input hashes and registration settings. The older run therefore cannot validate the accelerated kernel's whole-pipeline output.

The output checkpoint is under `build/plus-labeled-folded-numba-onepass-2026-10-12/`; `plus-full-workflow.json` is the restart manifest and `ribbon-quality.json` records the audit. The audit failed acceptance: GM-between-surfaces support **75.716%**, WM-inside-inner support **91.280%**, closest-surface thickness p05/median/p95 **0.170/1.394/5.325 mm**, and each mesh has 29 components (Euler characteristic 40). Level-set mean thickness is 4.247 mm and level-set nesting passes. All are synthetic results; no real-human PLUS validation is implied.

## Next

Keep the Python backend available for same-state checks. Before claiming whole-run equivalence or keeping the integrated speed result as a controlled benchmark, investigate why the registered FLAIR output changed between these runs, rerun both backends from the same persisted workflow checkpoint, and retain the required ribbon-quality gates. The integrated PLUS case remains nonterminal and fails tissue support and mesh quality.
