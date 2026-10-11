# P41z exact-crop source localization — 2026-10-11

## Evidence on the corrected iteration-6 surface checkpoint

This is a held-fixed selected-source comparison on the corrected 81×101×77
synthetic plus case, RAS crop 73×91×63. It does not compare the complete
MaCRUISEplus workflow with Java and does not validate real MRI. The Python code
under comparison was at commit `0929c81f`. The ACE/GVF and TGDM JSON reports
have `implementation_signature: null`, so the tested revision is pinned by this
report rather than by those machine-readable signatures.

| Stage | Java/Python result | Interpretation |
|---|---|---|
| TOADS topology correction and MaWM erosion/refinement | All 3 outputs bitwise exact over 418,509 voxels | No mismatch in these selected topology operations on this crop. Java RSS was below the watchdog measurement resolution. |
| ACE mode 3 | Enhanced GM, skeleton, and thinned skeleton bitwise exact over 418,509 voxels each; peak 0.250 GB under 8 GB | No mismatch in the selected ACE outputs. |
| GVF, 200 iterations | Mean absolute error 1.89e-8; max 1.86e-7; peak 0.309 GB under 8 GB | Within the pre-existing Java comparison limits (mean ≤1e-5, max ≤5e-5). |
| TGDM fields | All six fields have zero sign changes; five surface Dice scores are 1.0; surface mean distances are 7.58e-9–1.09e-7 mm and Hausdorff distances 7.62e-6–4.42e-5 mm | All agreed per-surface thresholds pass on the exact P41z crop. Values differ at some voxels only by small floating-point rounding. |
| Raw MaCRUISEplus boundary | Java equals Python, but both produce zero foreground voxels on the saved P41z crop | Exact but trivial for this case; not a useful nontrivial boundary validation. The companion folded-boundary fixture has incompatible dimensions/affine for this workflow. |

Machine-readable reports: [topology](p41z-topology-java-comparison-2026-10-11.json),
[ACE/GVF](p41z-ace-gvf-java-comparison-2026-10-11.json),
[TGDM](p41z-tgdm-java-comparison-2026-10-11.json), and the
[boundary diagnostic](p41z-boundary-java-diagnostic-2026-10-11.json).

## Where the ribbon failure appears

The quality audit reports that 67,674 of 67,677 input MaSeg cortical-GM voxels
are outside the final outer surface, and only 59.77% of the MaSeg WM voxels are
inside the inner surface. On the exact saved crop:

- The `multi_atlas_gm` map contains all 67,677 MaSeg label-100 voxels, and the
  `multi_atlas_wm` map contains all 37,699 label-44 voxels.
- All 67,677 label-100 voxels are positive in `initial_level_set`; none is in
  its nonpositive region. The initial `tissue_sign` map has value zero at all
  label-100 voxels.
- The raw boundary map is empty. The selected source only builds boundaries
  for labels ≥100 and suppresses boundaries touching a lower-valued label;
  this case has only one label ≥100 (label 100).

The positive initial level set at label-100 voxels is expected: MaCRUISEplus
seeds this field from topology-corrected TOADS WM, not discrete MaSeg cortical
labels. The empty one-parcel boundary is also expected because the source
suppresses candidate boundaries touching labels 1–99. Neither is evidence of
an implementation defect.

The zero fused `tissue_sign` is explained by a wiring mismatch upstream of
TGDM. Source MaCRUISE crops and passes `Atlas_BrainTissueMask` to
`CombineMAandTOADS`; this includes GM, WM, and other brain tissue. Python
`plus_full_workflow.py` currently passes `Atlas_Mask`, the WM-only mask. In
this saved crop all 67,677 label-100 voxels have zero fusion-mask support, even
though atlas GM is 255 at all of them and transferred TOADS GM is nonzero at
67,672. Consequently all 67,677 voxels take the outside-mask fusion branch;
none take confident-WM, confident-GM, or fallback branches, and fused GM/WM
are zero there. The selected topology, ACE, GVF, and TGDM method comparisons
above show no mismatch on their held-fixed inputs, but they cannot validate
the incorrect whole-workflow fusion input. The current run's surfaces remain
disconnected (2 components per
mesh, Euler −16) and too thin (median closest-point triangle gap 0.3725 mm;
level-set mean 0.8106 mm over 291 crossings). The 0.8106-mm statistic is sparse
and must not be treated as dense ribbon thickness.

The corrected registered loop was also continued from iteration 6 to 7. It
completed in **184.25 s** at **0.780 GB** peak sampled RSS, with energy ratio
**0.1319265813** and no stop. This is still far above the selected `1e-4`
threshold. It has not been composed through surfaces.

## Reproduction and scope

Generated intermediates are under the ignored
`build/p41z-case-normalization/` directory. ACE/GVF and TGDM inputs were
reconstructed from the saved Python crop and selected settings, serialized
once, and passed unchanged to their Java harnesses; Python products were
compared on that same grid and those same reconstructed inputs. This is a
held-fixed selected-method comparison, not an independently generated Java
workflow. The ACE/GVF and TGDM reports have a null implementation signature;
the tested source revision is pinned above.

- Recompare topology products: `.venv/bin/python scripts/compare_plus_topology_java.py --reference-dir build/p41z-case-normalization/topology-java --report /tmp/p41z-topology-recheck.json`.
- Recompare ACE/GVF products: `.venv/bin/python scripts/compare_plus_ace_gvf_java.py --reference-dir build/p41z-case-normalization/ace-gvf-java --report /tmp/p41z-ace-gvf-recheck.json`.
- Recompare TGDM products: `.venv/bin/python scripts/compare_plus_tgdm_java.py --workflow-root build/p41z-case-normalization/iteration-0006-surfaces/plus-surfaces --reference-dir build/p41z-case-normalization/tgdm-java-inputs --report /tmp/p41z-tgdm-recheck.json`.
- Regenerate Java products by following `scripts/compare_plus_topology_saved_java.sh`, `scripts/compare_plus_ace_gvf_saved_java.sh`, and `scripts/compare_plus_tgdm_saved_java.sh` with the saved crop and output directories. A nontrivial folded-boundary check is still unavailable on this crop because the committed folded fixture is not grid-compatible.

The latest full pytest suite remains **1,531 passed, 0 failed, 0 skipped** on
older implementation commit `9cade54`; this turn ran the focused lesion
regression group (**17 passed, 0 failed, 0 skipped in 0.86 s**) on the
corrected code. These new source comparisons are standalone harness results,
not a full pytest run.

## Next checkpoint

Use `Atlas_BrainTissueMask` instead of the WM-only `Atlas_Mask` in the Python
PLUS fusion input, matching the Java wrapper. First hold the current crop and
other fusion inputs fixed and compare Java/Python fusion products and branch
counts; then invalidate and resume the PLUS stages that depend on that input,
and recheck ribbon quality. Preserve the current case and thresholds as
diagnostics; do not relax them. A small two-parcel boundary fixture (labels 100
and 101, plus a paired lower-label suppression case) can separately validate
nontrivial source boundary behavior.
