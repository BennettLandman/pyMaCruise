# Human BASE attempt-02 TGDM held-fixed Java comparison — 2026-10-16

## Evidence category and scope

This is a selected-source, held-fixed **TGDM stage comparison** on the exact
human-derived Python crop from the regenerated attempt-02 BASE run. It is not a
fresh full MaCRUISE Java/MATLAB execution, nor evidence that upstream Python
inputs match those of the original run. The Java comparison consumes the same
five scalar arrays and three-component GVF produced by Python, then applies the
selected `MaGenericTGDM` nested-surface call sequence and V3 baseline settings.

## Inputs and source

- Python run: `build/sample-singularity-base-source-profile-2026-10-15/`
- Crop: 140×169×134 voxels; inclusive full-grid bounds `(21:160, 33:201,
  71:204)`; 1 mm axis-aligned RAS grid.
- Harness inputs, in Java loop order `(x,y,z)`: `wmfill_topology_corrected`,
  `wmfill_cropped`, `gm_cropped`, `ace_enhanced_gm`,
  `ace_thinned_skeleton`, and the three-component `gvf`.
- The comparison Java source is
  [`PlusTGDMReferenceHarness.java`](../../tests/java_reference/PlusTGDMReferenceHarness.java),
  which calls the selected original
  [`MaGenericTGDM.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/plugins/CRUISE/tgdm/MaGenericTGDM.java)
  implementation. Although this harness was first added for PLUS work, the
  sequence and parameters used here are the BASE `tgdm_nested_images` sequence;
  no PLUS-specific registration or preprocessing is run.
- Source hashes: harness `a5b7c53b94d0d70a68999ec8c35b8985a5c1a8b2b114fe12be07afee884a4a3a`;
  MaGenericTGDM `66124d4180864edc77f6fff280ab817eb0fffc061fdf427c9d8d725f850b211d`;
  selected distance-field source `a0df8859367a52a6273c488fed74d3fbdf0bbf5776b6e6cfea15261a5b290811`.
- Selected settings: initial iso 126.9, WM iso 0.5, GM iso 0.4; initial
  smoothing 6 iterations; inner 4, central 7, central shrink 1, and outer 7;
  connectivity 18/6. Java uses `-Xmx6g` inside a 16 GB RSS watchdog.

## Result

Python's crop outputs were compared to Java's central, inner, and outer
zero-level surfaces in physical RAS millimeters. All agreed surface tolerances
pass: Dice ≥0.97, symmetric mean nearest-surface distance ≤0.01 mm, and
sampled Hausdorff distance ≤1 mm.

| Surface | Dice at zero | Symmetric mean distance (mm) | Hausdorff (mm) | Sign-different voxels |
|---|---:|---:|---:|---:|
| Central | 1.000000 | 0.0000000035 | 0.0000152 | 0 |
| Inner | 1.000000 | 0.0000000042 | 0.0000149 | 0 |
| Outer | 0.999972 | 0.0003529 | 0.710120 | 55 |

The inner and central Java/Python meshes have identical vertex and face counts.
The outer meshes differ by six vertices and twelve faces, while still passing
the agreed position tolerances. Full level-set values are not bit-identical:
104,423 central, 107,195 inner, and 2,160,842 outer values differ. Mean
absolute level-set errors are 0.03394, 0.02463, and 0.05474, respectively;
maximum errors are 14.51, 13.42, and 15.59. Therefore this report supports
surface-position agreement, not full-field numeric identity.

The Java replay completed successfully. Peak sampled RSS was **1.317 GB**;
the 16 GB watchdog did not trigger. Raw arrays, Java outputs, and machine
metrics remain in the ignored `build/.../java-tgdm-attempt-02/` directory.
The compact committed [machine-readable comparison](in-vivo-root-sample-base-attempt-02-tgdm-2026-10-16.json)
records source/input/output SHA-256 values, parameters, resource use, and the
full per-surface metrics.

## Interpretation

This closes one additional selected-source comparison on a real-human-derived
BASE crop. The numerical target is the Python-generated handoff and the chosen
surface-position criteria. It does not close whole-workflow equivalence,
confirm the unknown provenance/settings of the downloaded historical result,
or replace a matched original MATLAB/MIPAV run.
