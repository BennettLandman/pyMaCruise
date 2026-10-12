# P41z iteration-15 selected TGDM source comparison — 2026-10-11

## Result

The selected Java `MaGenericTGDM` wrapper sequence passes all five agreed
surface gates when compared with the Python TGDM outputs from the iteration-15
P41z held-fixed workflow. The Java harness ran under an 8 GB RSS watchdog and
peaked at **0.626 GB**. The input crop is 73×91×63 on an approximately 2 mm
RAS grid.

| Field | Sign-different voxels | Dice at zero | Symmetric mean distance (mm) | Hausdorff (mm) | Gates |
| --- | ---: | ---: | ---: | ---: | --- |
| Smoothed initial | 0 | 1.000000 | 8.86e-8 | 1.58e-5 | Pass |
| Initial inner | 0 | 1.000000 | 9.04e-9 | 7.66e-6 | Pass |
| Central | 0 | 1.000000 | 8.19e-9 | 7.64e-6 | Pass |
| Inner | 0 | 1.000000 | 8.81e-9 | 7.64e-6 | Pass |
| Outer | 7 | 0.999970 | 0.001688 | 0.494472 | Pass |

Four of the six saved TGDM products (including the initial level set) have no
sign changes. Small floating-point differences remain, and the outer field has
seven changed signs; every surface remains within the 0.97 Dice, 0.01 mm
symmetric mean-distance, and 1 mm symmetric Hausdorff limits.

## Evidence scope and interpretation

This is a synthetic selected-stage comparison, not full Java workflow
equivalence. The Java harness reuses the hash-checked Python signed-distance
seed; tissue, ACE, and GVF inputs are reconstructed on the same crop grid; and
TOADS/LesionTOADS remain held fixed. This result does not fix the iteration-15
anatomy-quality failures: GM support is 92.890%, mesh Euler characteristic is
−18, and production surface-thickness percentiles remain out of bounds. The
registered loop is also nonterminal at energy ratio 0.13346674. See the
[iteration-15 ribbon audit](p41z-iteration15-scaled-wrapper-ribbon-audit-2026-10-11.md)
and [loop continuation report](p41z-registered-loop-iterations-12-to-15-2026-10-11.md).

## Reproduction

```sh
scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/iterations-0012-to-0015-surfaces/plus-surfaces \
  build/p41z-case-normalization/iterations-0012-to-0015-surfaces/tgdm-java \
  build/p41z-case-normalization/iterations-0012-to-0015-surfaces/tgdm-java-comparison.json \
  build/p41z-case-normalization/iterations-0012-to-0015-surfaces/plus-surfaces/run.json
```

The machine-readable result and raw Java watchdog log are in the ignored
`build/p41z-case-normalization/iterations-0012-to-0015-surfaces/` directory.
