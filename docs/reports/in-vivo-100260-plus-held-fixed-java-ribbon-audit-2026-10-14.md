# Human PLUS held-fixed Java ribbon audit: case 100260

**Evidence category:** real-human anatomy diagnostics comparing the selected Java TGDM sequence with Python on the same crop-ready workflow inputs. This is not an independently executed Java/MATLAB PLUS workflow.

## Result

The selected Java TGDM replay reproduces the current Python ribbon-quality concerns on the same case and inputs. Both produce the same 121 connected, watertight surface components with Euler characteristic 70; WM-inside support is 97.63%, GM-between support is 43.49%, and median closest-point thickness is about 1.063 mm. The Java and Python TGDM fields have identical signs except for 10 outer-field voxels, and the resulting ribbon metrics differ only slightly.

| Metric | Selected Java TGDM | Python TGDM |
|---|---:|---:|
| WM inside inner surface | 97.6336% | 97.6336% |
| GM between inner and outer surfaces | 43.4880% | 43.4883% |
| Inner outside central voxels | 46 | 46 |
| Central outside outer voxels | 0 | 0 |
| Surface connected components (inner / central / outer) | 121 / 121 / 121 | 121 / 121 / 121 |
| Euler characteristic (inner / central / outer) | 70 / 70 / 70 | 70 / 70 / 70 |
| Closest-point median thickness | 1.062714 mm | 1.062731 mm |
| Closest-point p05 / p95 | 0.243392 / 2.528595 mm | 0.243379 / 2.528643 mm |
| Level-set mean thickness | 2.206767 mm | 2.206755 mm |

This indicates that the low GM ribbon support, thin ribbon, and disconnected topology are already present in the selected Java TGDM result when given the Python workflow's identical crop inputs. It therefore does not identify a Python-vs-Java TGDM defect as the cause. It also does not determine whether the source workflow's independent upstream products (registration, LesionTOADS, Dura transfer, crop, or topology inputs) differ. The broader source-workflow comparison remains open.

The comparison used the hash-verified cropped MaSeg input and the Java fields saved by the corrected held-fixed comparison. All inputs and generated fields were on a shared RAS grid of `140×178×134`. The audit peaked at **1.241 GB** under an **8 GB** RSS watchdog. Machine-readable metrics, including per-field differing voxel counts and the full surface summaries, are in [the JSON report](in-vivo-100260-plus-held-fixed-java-ribbon-audit-2026-10-14.json).

## Reproduction

From the repository root, after generating the Java held-fixed replay described in the [comparison report](in-vivo-100260-plus-held-fixed-java-comparison-2026-10-14.md):

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/audit_plus_held_fixed_java_ribbon.py \
  --reference-dir build/in-vivo-100260-plus-java-corrected-replay-2026-10-14 \
  --provenance build/in-vivo-100260-plus-java-corrected-replay-2026-10-14/input-provenance.json \
  --report docs/reports/in-vivo-100260-plus-held-fixed-java-ribbon-audit-2026-10-14.json
```

This is an anatomy diagnostic only. The Java replay holds the Python-generated pre-TGDM workflow inputs fixed and must not be described as an independent or complete PLUS source comparison.
