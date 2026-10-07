# Folded-cortex TGDM comparison against Java — 2026-10-07

## Scope

This is a synthetic, stage-wise source comparison of the default Python TGDM
pipeline with the supplied `MaGenericTGDM.java` CVS 1.5 and
`DistanceFieldyh.java` CVS 1.1. It uses the committed 41×41×41 analytic folded
cortex inputs in `tests/data/synthetic/folded_cortex/`; all five input arrays
are read from those exact NIfTIs by both implementations. This does not claim
whole MaCRUISE workflow equivalence or real-MRI validation.

The Java reference harness now accepts an input directory in addition to its
original no-argument synthetic sphere mode. It supports the fixture's 3D
little-endian uint8 and float32 NIfTIs, maps their x/y/z voxels to the JIST
array order, and writes RAS identity-affine stage outputs. The harness runs the
same selected stage parameters as the existing 25³ reference. Java ran with a
2 GiB heap cap.

## Reproduction

```sh
scripts/generate_java_tgdm_folded_reference.sh
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --report-interval-seconds 15 -- \
  .venv/bin/python scripts/compare_tgdm_java_reference.py \
  --input-root tests/data/synthetic/folded_cortex \
  --reference-root tests/data/synthetic/folded_cortex/java_reference \
  --report docs/reports/tgdm-folded-cortex-java-comparison-2026-10-07.json
```

The Java harness completed all stages and wrote checksummed `.nii.gz` references. The Python comparison completed with
exit 0 on committed implementation `89e7c71`; sampled peak RSS was **0.168 GB** under the 8 GB watchdog. The harness's original 25³ sphere mode was also rerun to a temporary directory and passed
the existing Java comparison thresholds, preserving backward compatibility. Inputs and
outputs are RAS on a 1 mm identity grid. The tissue-sign field was bit-exact.
Every evolved level set had mask Dice 1.0, symmetric mean surface distance at
or below 0.003108 mm, and symmetric Hausdorff distance at or below 0.259171 mm.
All are within the project acceptance thresholds (Dice ≥0.97, mean distance
≤0.01 mm, Hausdorff ≤1 mm).

| Stage | Dice | Symmetric mean surface distance (mm) | Hausdorff (mm) | Voxel MAE |
|---|---:|---:|---:|---:|
| Signed seed | 1.000 | 0.000000 | 0.000000 | 3.60e-8 |
| Smoothed initial | 1.000 | 0.000099 | 0.001960 | 0.000181 |
| Initial inner | 1.000 | 0.000001 | 0.000005 | 0.000001 |
| Central | 1.000 | 0.000052 | 0.000824 | 0.001142 |
| Central shrink | 1.000 | 0.000387 | 0.006542 | 0.000312 |
| Final inner | 1.000 | 0.0000005 | 0.000002 | 0.000001 |
| Outer | 1.000 | 0.002852 | 0.259171 | 0.003731 |

See the machine-readable [full metric report](tgdm-folded-cortex-java-comparison-2026-10-07.json).

## Interpretation

This extends the existing 25³ concentric-sphere comparison to a larger folded
surface with narrow-band geometry. It supplies direct Java comparison evidence
for all nested TGDM level-set stages on this fixture. The outer surface has the
largest numerical deviation but remains within the agreed surface thresholds.
The 25³ comparison had Dice 1.0 at every level-set stage, maximum symmetric
mean surface distance 0.000652 mm and Hausdorff 0.018173 mm. The result supports
TGDM stage behavior on these two synthetic geometries; it
does not establish convergence or numerical equivalence for arbitrary images,
other parameter profiles, GVF inputs, wrapper composition, or the complete
MaCRUISE run.
