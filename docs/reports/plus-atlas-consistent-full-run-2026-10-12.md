# Atlas-aligned synthetic PLUS full run — 2026-10-12

## Run status and validation

The committed-input candidate in `tests/data/synthetic/plus_atlas_consistent/` was run through the public raw-input `plus-full-run` command with source-faithful zero-variance handling and `--max-passes 1`. All nine public workflow stages completed and wrote restartable manifests, the Dura outputs, native-grid transfers, three RAS-mm VTK surfaces, three RAS level sets, and thickness products. The LesionTOADS loop stopped at the requested cap at iteration 2 (energy ratio 0.108983), so the result is a complete one-pass checkpoint, not a terminal workflow.

An independent manifest walk verified **95/95** reachable output/input hashes, **69/69 finite RAS NIfTIs**, and no hash or orientation errors. Peak sampled RSS was **5.725 GB** under a 16 GB watchdog. This is synthetic integration evidence only; the case is derived from the bundled atlas and is not independent truth for atlas-registration or TOADS accuracy.

## Ribbon audit

The explicit audit fails five quality gates. MaSeg-like support is WM-inside **99.989%**, GM-between **66.709%**, and CSF-outside **99.648%** (required minimum 97%). The inner surface lies outside the central surface at two voxels. Each surface mesh is watertight and manifold but has **17 connected components** and Euler characteristic **−22**. Production closest-surface thickness has mean/median/p05/p95 **4.059/4.853/0.161/7.065 mm**; the very thin p05 fails the 2 mm floor. Level-set subtraction mean thickness is **7.075 mm**, beyond the 5 ± 1 mm bound. The loop remains nonterminal.

This atlas-derived fixture produces different support behavior from the earlier 65³ folded fixture: WM support is higher, while GM support is lower. Label/atlas alignment alone does not yield an acceptable ribbon, and this fixture cannot serve as independent validation because its inputs were built from the same atlas used by the workflow. Keep both fixtures and audits as regression evidence, and develop an independent folded phantom with planted known surfaces and explicit lesion/no-lesion selection before treating synthetic ribbon acceptance as established.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/python tests/data/synthetic/plus_atlas_consistent/generate_fixture.py
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_atlas_consistent/t1.nii.gz \
  --flair tests/data/synthetic/plus_atlas_consistent/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_atlas_consistent/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_atlas_consistent/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-atlas-consistent-full-baseline-2026-10-12 \
  --max-passes 1
PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-atlas-consistent-full-baseline-2026-10-12/plus-full-workflow.json \
  --report build/plus-atlas-consistent-full-baseline-2026-10-12/ribbon-quality.json \
  --allow-checkpoint
```

The final run manifest and full audit are local `build/` outputs. The compact committed audit is [plus-atlas-consistent-ribbon-quality-2026-10-12.json](plus-atlas-consistent-ribbon-quality-2026-10-12.json); the generated case and assumptions are documented in the fixture directory.
