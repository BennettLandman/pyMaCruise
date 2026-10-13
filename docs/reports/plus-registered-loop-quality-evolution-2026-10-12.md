# Registered synthetic PLUS loop quality evolution — 2026-10-12

## Scope and evidence category

This report tracks continued LesionTOADS registered-loop iterations on the committed `plus_atlas_consistent` synthetic input. It is **synthetic workflow evidence**, not proof of Java equivalence or real-human accuracy. The fixture's tissue labels were grouped from the same atlas used by the workflow, so it is not independent truth. Every checkpoint below is nonterminal; surface composition is explicitly marked as checkpoint composition.

The iteration-2→3 replay also compared the full atlas-sized registered shape-prior output to the selected Java harness. All **71,127,680 float32 values matched exactly**, as did selected thinning, statistics, membership, registration, growth, relation, and saved transition outputs. This validates that source-selected operation on this case, not the whole MaCRUISEplus workflow.

## Progress and quality

| Iteration | Energy ratio | WM inside inner | GM between inner/outer | CSF outside outer | Surface components / Euler | Surface median / p05 / p95 (mm) | Level-set mean (mm) | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 2 | 0.108983 | 99.989% | 66.709% | 99.648% | 17 / −22 | 4.853 / 0.161 / 7.065 | 7.075 | Fails ribbon gates |
| 7 | 0.002044 | 99.954% | 91.505% | 99.630% | 16 / −46 | 3.956 / 0.167 / 8.514 | 8.571 | Fails ribbon gates |
| 10 | 0.000511 | 99.955% | 91.305% | 99.631% | 23 / −28 | 3.958 / 0.167 / 8.508 | 8.553 | Fails ribbon gates |
| 12 | 0.000956 | 99.955% | 91.352% | 99.630% | 21 / −30 | 3.959 / 0.167 / 8.512 | 8.558 | Fails ribbon gates |

The required label-support minimum is 97%; the surface p05 minimum is 2 mm, median must be within 1 mm of 5 mm, p95 must be at most 8 mm, level-set mean must be within 1 mm of 5 mm, and each mesh must be one closed connected genus-zero component. Every measured checkpoint fails multiple gates. WM and CSF support pass, but GM support remains near 91.3% after iteration 7. The energy ratio is not monotonic: it drops through iteration 11 to 0.000513, then rises to 0.000956 at iteration 12. Small energy change alone is not an acceptance criterion.

The iteration-12 meshes are watertight and manifold, but each has 21 connected components and Euler characteristic −30. Inner and central voxel nesting disagrees at four voxels. Thickness tails are also persistent: closest-surface p05 is about 0.167 mm and p95 about 8.512 mm; level-set mean is 8.558 mm versus the 5 mm target.

## Runtime and memory

- Full raw-input synthetic `plus-full-run`, one requested loop pass: 5.725 GB peak sampled RSS; all nine stages completed, 95/95 reachable manifest hashes verified, and 69/69 NIfTIs were finite and RAS.
- Selected Java iteration-2→3 comparison: 9.027 GB peak sampled RSS under a 16 GB watchdog.
- Continuation through iterations 5–7: 8.227 GB peak sampled RSS.
- Continuation through iterations 8–10: 7.280 GB peak sampled RSS.
- Continuation through iterations 11–12: 10.440 GB peak sampled RSS.
- Iteration-10/12 surface composition: 3.444/3.449 GB peak sampled RSS; audits stayed below 0.8 GB.

All continuations completed without crossing the 16 GB watchdog limit. The Numba registered-shape-prior path reduced the large-case sampling phase from about 224 seconds to under one second, while the whole loop remains dominated by source-selected optimization and growth stages.

## Reproduction

The raw-input baseline command and fixture provenance are in [the atlas-aligned full-run report](plus-atlas-consistent-full-run-2026-10-12.md). Continue from the iteration-10 checkpoint with:

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/plus-atlas-consistent-loop-iterations8to10-2026-10-12/lesion-registered-loop.json \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-atlas-consistent-loop-iterations11to12-2026-10-12 \
  --iterations 99 --thinning-limit 1000 --max-passes 2
```

Compose the nonterminal checkpoint with `scripts/compose_terminal_synthetic_plus.py --allow-checkpoint` and audit it with `scripts/audit_plus_ribbon_quality.py --allow-checkpoint`. The detailed machine-readable iteration-10 and iteration-12 audits are retained in ignored `build/` output for local reproduction; compact iteration-12 findings are recorded above.

## Interpretation and next milestone

This case shows a substantial increase in GM support from iteration 2 to iteration 7, followed by a plateau through iteration 12 while surface topology and thickness tails remain poor. Continuing the same iterations alone is unlikely to establish an acceptable synthetic ribbon. The next implementation work should diagnose the GM-label/surface disagreement and topology/thickness failures against intermediate stage outputs, then create an independent planted-label folded phantom so quality gates can distinguish algorithm defects from atlas-derived fixture bias. Do not describe the current PLUS workflow as fully validated until those checks pass and a terminal source-selected loop is demonstrated.
