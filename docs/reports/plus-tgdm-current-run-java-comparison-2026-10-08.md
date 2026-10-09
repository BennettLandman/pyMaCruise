# Held-fixed MaCRUISEplus TGDM Java comparison — 2026-10-08

## Scope

This is a selected-source comparison of the nested MaCRUISEplus TGDM stage on the exact held-fixed inputs reconstructed from the saved 64×80×64 folded synthetic plus workflow. The crop-ready TGDM grid is **66×76×66**. Upstream topology correction, ACE, GVF, TOADS, LesionTOADS, registration, and crop inputs are held fixed; this does not claim whole-workflow Java/MATLAB equivalence or real-MRI validation.

The Python comparison utility reruns the saved plus TGDM stage from the stage's exact crop inputs and settings before launching Java. It first requires all six Python outputs (including the tissue-sign image) to match the saved stage arrays exactly. They did. The selected Java harness then runs the `MaGenericTGDM` wrapper sequence on the same initial-membership, WM, GM, ACE-GM, thin-skeleton, and GVF arrays.

## Reproduction

From the repository root:

```sh
./scripts/compare_plus_tgdm_saved_java.sh
```

The default workflow is `build/plus-pyramid-source-current/natural-run/plus-surface`. Input hashes, workflow manifest hash/signature, and the held-fixed binary hash are recorded in the generated ignored build artifact `build/plus-tgdm-current-java-2026-10-08/input-provenance.json`. The result metrics are saved in [the machine-readable comparison](plus-tgdm-current-run-java-comparison-2026-10-08.json).

Java exited successfully under an **8 GB RSS watchdog**, with **0.442 GB peak sampled RSS** and no RSS termination. The maximum heap was 6 GB. The exact held-fixed arrays were replayed before comparison; five Java/Python scalar output arrays differ in floating-point values, but no sign-classification voxels differ.

## Results

| Output | Changed voxels | Maximum absolute voxel error | Sign changes | Dice at zero | Symmetric mean surface distance | Hausdorff |
|---|---:|---:|---:|---:|---:|---:|
| Initial level set | 9,220 | 9.54×10⁻⁷ | 0 | — | — | — |
| Initial inner | 56,918 | 5.73×10⁻⁶ | 0 | — | — | — |
| Central | 59,495 | 0.00442 | 0 | 1.0 | 1.65×10⁻⁶ mm | 0.000347 mm |
| Inner | 61,133 | 0.000215 | 0 | 1.0 | 2.19×10⁻⁷ mm | 0.0000866 mm |
| Outer | 59,154 | 0.0485 | 0 | 1.0 | 1.69×10⁻⁶ mm | 0.000328 mm |

All three delivered surfaces meet the agreed acceptance thresholds: Dice ≥0.97, symmetric mean nearest-surface distance ≤0.01 mm, and Hausdorff distance ≤1 mm. The outer scalar field's maximum voxel error is localized in values without changing the zero-level sign mask; the geometric output remains within tolerance.

## Interpretation

This is synthetic, held-fixed Java source evidence for the plus TGDM sequence, not an end-to-end workflow comparison. It narrows the P37 gap for the folded plus case. The upstream MaSeg/topology/ACE/GVF sequence, LesionTOADS registration and convergence, cross-grid sensitivity, and the complete plus composite still need source-aligned evidence. Human-MRI validation also remains open and requires a matching user-provided brain mask.
