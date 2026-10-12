# P41z TGDM Java comparison after Dura wrapper-scale correction — 2026-10-11

## Result

The selected Java `MaGenericTGDM` wrapper sequence passes the agreed surface
gates on the corrected P41z held-fixed inputs. The Python implementation
replayed all six saved TGDM products exactly before generating the Java input
arrays. Java and Python surfaces have identical zero-level Dice for the first
four fields and differ at only three sign voxels in the outer field; all five
surface comparisons pass the Dice ≥0.97, symmetric mean distance ≤0.01 mm,
and symmetric Hausdorff ≤1 mm limits.

| TGDM field | Sign-different voxels | Dice at zero | Symmetric mean distance (mm) | Hausdorff (mm) | Gates |
| --- | ---: | ---: | ---: | ---: | --- |
| Smoothed initial | 0 | 1.000000 | 8.89e-8 | 2.19e-5 | Pass |
| Initial inner | 0 | 1.000000 | 9.85e-9 | 1.14e-5 | Pass |
| Central | 0 | 1.000000 | 6.42e-9 | 8.87e-6 | Pass |
| Inner | 0 | 1.000000 | 8.63e-9 | 7.64e-6 | Pass |
| Outer | 3 | 0.999987 | 0.001531 | 0.444017 | Pass |

The full comparison JSON is retained in the ignored build output at
`build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/tgdm-java-comparison.json`.
Java peak sampled RSS was **0.622 GB** under the 8 GB watchdog. The RAS crop
was 73×91×63 at approximately 2.01 mm spacing.

## Interpretation and limits

The corrected Python surfaces fail the separate ribbon-quality acceptance
despite passing Java surface equivalence: GM label-100 support between inner
and outer is **92.835%**, below the 97% criterion. The held-fixed Java TGDM
replay is subject to the same poor tissue support. This localizes the remaining
problem upstream of numerical TGDM parity for these particular inputs, but it
does not identify whether the cause is the held-fixed TOADS/LesionTOADS state,
nonterminal iteration-11 loop, seed construction, or their interaction.

The Java replay uses the saved Python initial signed-distance field and
reconstructed same-grid tissue, ACE, and GVF inputs. TOADS and LesionTOADS are
held fixed; this is a selected-stage comparison, not a full Java MaCRUISEplus
run, whole-loop equivalence, or real-human-MRI acceptance.

## Reproduction

```sh
scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/plus-surfaces-corrected \
  build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/tgdm-java \
  build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/tgdm-java-comparison.json \
  build/p41z-case-normalization/iterations-0008-to-0011-scaled-wrapper/plus-surfaces-corrected/run.json
```
