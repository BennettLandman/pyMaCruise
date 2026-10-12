# P41z iteration-19 selected TGDM Java comparison — 2026-10-11

## Result

The selected Java `MaGenericTGDM` wrapper sequence passes all five agreed
surface gates when compared with the iteration-19 Python plus workflow. Four
surfaces have zero sign-changing voxels. The outer surface has six changed
signs but remains within all thresholds. The Java harness peaked at **0.622
GB** under an 8 GB watchdog.

| Field | Sign-different voxels | Dice at zero | Symmetric mean distance (mm) | Hausdorff (mm) | Gates |
| --- | ---: | ---: | ---: | ---: | --- |
| Smoothed initial | 0 | 1.000000 | 8.46e-7 | 0.001964 | Pass |
| Initial inner | 0 | 1.000000 | 9.63e-8 | 0.000709 | Pass |
| Central | 0 | 1.000000 | 1.42e-8 | 1.53e-5 | Pass |
| Inner | 0 | 1.000000 | 7.85e-9 | 7.67e-6 | Pass |
| Outer | 6 | 0.999975 | 0.001899 | 0.421216 | Pass |

The six saved Python products used to prepare the Java inputs replayed exactly
in Python. Floating-point differences remain between Java and Python arrays,
but all compared surfaces meet Dice ≥0.97, symmetric mean distance ≤0.01 mm,
and symmetric Hausdorff ≤1 mm.

## Evidence scope and interpretation

This is a selected-stage synthetic comparison, not full Java registered-loop
or whole-workflow equivalence. The Java harness reuses the saved Python
signed-distance seed, reconstructs tissue/ACE/GVF inputs on the same crop
grid, and holds TOADS/LesionTOADS fixed.

The comparison does not resolve the iteration-19 ribbon-quality failure:
GM-between-surfaces support is 92.906%, each mesh Euler characteristic is
−16, and production surface-thickness median/p05/p95 is 6.505/0.236/12.321
mm. The registered loop is still nonterminal at ratio 0.13316543. See the
[iteration-16-to-19 checkpoint and ribbon audit](p41z-iteration16-to-19-checkpoint-audit-2026-10-11.md).

## Reproduction

```sh
scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/plus-surfaces \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/tgdm-java \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/tgdm-java-comparison.json \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/plus-surfaces/run.json
```

The machine-readable comparison and Java watchdog log are in the ignored
`build/p41z-case-normalization/iterations-0016-to-0019-surfaces/` directory.
