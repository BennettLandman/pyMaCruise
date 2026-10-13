# Correctly labeled atlas-aligned synthetic PLUS run — 2026-10-12

## Evidence category and source label mapping

This is a **synthetic integration and source-selected comparison**, not real-human validation. The fixture uses the bundled topology atlas for anatomy and is not independent truth for atlas registration or TOADS.

A source-derived MaSeg grouping correction was required before interpreting the ribbon audit: topology class 15 alone maps to cortical GM 100; classes 25–27 map to cortical WM 44; classes 5–6 map to CSF 5; topology 1 maps to background 0; classes 14, 16–18, 22, and 24 map to other tissue 75. A regression test locks this mapping. The earlier atlas-aligned reports used an incorrect mapping and are superseded; see the erratum in [the earlier report](plus-atlas-consistent-full-run-2026-10-12.md) and [the registered-loop history](plus-registered-loop-quality-evolution-2026-10-12.md).

## Runnable workflow and completion

The corrected inputs were processed with the public `plus-full-run` CLI through all nine stages, initially capped at one loop pass. The registered loop was then resumed from the saved checkpoint and stopped naturally by convergence at iteration 10 with energy ratio **0.0000464068**. Terminal surfaces were composed from that saved state and explicitly marked as the terminal composition.

```sh
PYTHONPATH=src .venv/bin/python tests/data/synthetic/plus_atlas_consistent/generate_fixture.py
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_atlas_consistent/t1.nii.gz \
  --flair tests/data/synthetic/plus_atlas_consistent/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_atlas_consistent/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_atlas_consistent/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-atlas-consistent-corrected-full-2026-10-12 \
  --max-passes 1
```

Resume the loop from the stage-1 manifest, then from iteration 7. The second continuation stops at iteration 10 before its five-pass cap because the source-selected convergence test passes.

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest \
    build/plus-atlas-consistent-corrected-full-2026-10-12/lesion_workflow/registered-loop/lesion-registered-loop.json \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-atlas-consistent-corrected-loop-iterations3to7-2026-10-12 \
  --iterations 99 --thinning-limit 1000 --max-passes 5

env PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest \
    build/plus-atlas-consistent-corrected-loop-iterations3to7-2026-10-12/lesion-registered-loop.json \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-atlas-consistent-corrected-loop-iterations8to12-2026-10-12 \
  --iterations 99 --thinning-limit 1000 --max-passes 5

Compose the terminal surface products from the second continuation:

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/plus-atlas-consistent-corrected-loop-iterations8to12-2026-10-12/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_atlas_consistent/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_atlas_consistent/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_atlas_consistent/maseg.nii.gz \
  --output-dir build/plus-atlas-consistent-corrected-iteration10-compose-2026-10-12
```

The complete workflow manifest reports all nine stages complete. A recursive manifest audit verified **95/95** input/output hashes and **69/69** finite RAS NIfTIs. The corrected raw-input run peaked at **6.413 GB** sampled RSS; iterations 3–7 peaked at **10.214 GB**, iterations 8–10 at **10.018 GB**, terminal composition at **3.168 GB**, and the ribbon audit below **0.7 GB**, all under their watchdog limits.

The saved corrected terminal run passes held-fixed Java comparisons for all six TGDM level-set outputs and all four Dura-output channels. Every evolved TGDM surface passes zero-level Dice ≥0.97, mean surface distance ≤0.01 mm, and Hausdorff ≤1 mm; the outer surface is limiting (Dice 0.999943, mean 0.001837 mm, Hausdorff 0.937940 mm). The Dura output adapter matches exactly after the source's scale-by-255 operation. The historical JIST `UNCHANGED` Dura transfer also matches Python bit-for-bit, but preserves the cropped `(166, 206, 160)` extents while assigning the T1 reference header; the recommended physical-affine mode resamples to the `(172, 216, 164)` T1 grid. Upstream membership evolution and transfer-mode effects on the composed ribbon remain open. See [the held-fixed TGDM comparison](plus-atlas-consistent-corrected-tgdm-java-2026-10-12.md), [Dura output comparison](plus-atlas-consistent-corrected-dura-java-2026-10-12.md), and [JIST transfer comparison](plus-atlas-consistent-corrected-jist-transfer-java-2026-10-12.md).

## Ribbon acceptance

The explicit terminal audit still fails five quality gates:

| Check | Result | Required |
|---|---:|---:|
| WM inside inner surface | 99.955% | ≥97% |
| GM between inner and outer surfaces | 96.468% | ≥97% |
| CSF outside outer surface | 99.594% | ≥97% |
| Inner outside central nesting | 1 voxel | 0 |
| Surface mesh topology | 4 components, Euler −8 (each mesh) | 1 component, Euler 2 |
| Closest-surface thickness median / p05 / p95 | 4.200 / 0.150 / 8.479 mm | 4–6 / ≥2 / ≤8 mm |
| Level-set subtraction mean thickness | 8.623 mm | 4–6 mm |

GM support is only 0.532 percentage points below the agreed threshold, while the thickness tail and topology requirements remain material failures. Natural loop convergence does not imply ribbon acceptance.

## JIST surface extraction comparison

The selected JIST `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)` implementation and the production Python source-derived case-template extractor ran on the identical saved terminal level sets. For all three surfaces, the Java harness confirmed the exact level-set hash. Java JIST and the Python production extractor generated identical vertex/face counts, component count (4), Euler characteristic (−8), and unordered triangle connectivity (100% match). Symmetric sampled vertex-to-triangle mean distance was below **0.00000004 mm**, with sampled maximum at most **0.000245 mm**. Applying the same Python thickness operator gives mean thicknesses of 4.14115814 mm and 4.14115797 mm (delta **−0.00000016 mm**).

Thus, on this corrected terminal case, the four-component topology and ribbon-thickness failures originate in the level sets, not in mesh extraction: Java JIST and the Python production port produce the same mesh from the same fields. This does not settle triangulator behavior on other level sets; earlier source comparisons show that topology can differ on other fixtures. See [the exact-checkpoint JIST comparison](plus-terminal-jist-surface-comparison-2026-10-12.md) and its adjacent JSON. The production extractor uses the Java-derived 18/6 case-template table; the scikit-image comparison remains historical.

## Reproduction of the quality and surface checks

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/audit_plus_ribbon_quality.py \
  build/plus-atlas-consistent-corrected-iteration10-compose-2026-10-12/synthetic-plus-terminal-surface-composition.json \
  --report build/plus-atlas-consistent-corrected-iteration10-compose-2026-10-12/ribbon-quality.json

PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/compare_jist_surface_checkpoint.py \
  build/plus-atlas-consistent-corrected-iteration10-compose-2026-10-12/synthetic-plus-terminal-surface-composition.json \
  --output docs/reports/plus-terminal-jist-surface-comparison-2026-10-12.json \
  --markdown docs/reports/plus-terminal-jist-surface-comparison-2026-10-12.md
```

The next PLUS milestone is to diagnose the remaining GM boundary and thickness failures upstream of TGDM, using the corrected MaSeg labels, held-fixed TGDM comparison, and terminal JIST comparison as fixed references. Keep the independent planted-fold synthetic fixture as a separate acceptance case.
