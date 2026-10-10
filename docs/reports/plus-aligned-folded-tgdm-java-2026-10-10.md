# Atlas-aligned synthetic plus TGDM comparison — 2026-10-10

## Scope

This compares the selected `MaGenericTGDM` plus wrapper with the saved
Python TGDM stage from the terminal P41z composition. It uses the same
60×75×52 RAS crop, the hash-checked Python signed-distance seed, tissue fields,
ACE-GM/skeleton, and GVF. Tissue, ACE, and GVF input arrays were reconstructed
from the hash-checked crop inputs using the committed Python stages. The Java
process was run under an 8 GB RSS watchdog and peaked at **0.409 GB**.

This is a held-fixed TGDM source comparison only. It does not establish
whole-workflow equivalence, LesionTOADS convergence, or real-MRI behavior.

## Results

| Surface | Dice at zero | Symmetric mean nearest-surface distance | Hausdorff |
|---|---:|---:|---:|
| Central | 1.0 | 0.0000754 mm | 0.1008 mm |
| Inner | 1.0 | 0.0000211 mm | 0.0463 mm |
| Outer | 1.0 | 0.0000764 mm | 0.1008 mm |

All surfaces pass the agreed minimum Dice 0.97, maximum symmetric mean
distance 0.01 mm, and maximum Hausdorff 1 mm criteria. The scalar outputs
differ in floating-point values, but none of the zero-level sign masks differ.
This explains why the terminal ribbon audit's submillimeter thickness is not a
Python-versus-Java TGDM discrepancy for these held-fixed inputs. The audit's
label-support and thickness failures remain real; the next work should inspect
how upstream tissue/atlas fields and synthetic anatomy seed the selected
evolution.

## Reproduction

The raw arrays and Java class files are generated under ignored `build/` paths.
The committed machine-readable comparison is
[plus-aligned-folded-tgdm-java-2026-10-10.json](plus-aligned-folded-tgdm-java-2026-10-10.json).
The preparation path reuses hash-checked saved TGDM outputs and reconstructs
only the held-fixed input fields, avoiding a second multi-minute nested TGDM
run. Java compilation and execution completed successfully.
