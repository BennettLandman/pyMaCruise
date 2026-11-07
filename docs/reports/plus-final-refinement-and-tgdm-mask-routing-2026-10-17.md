# PLUS final refinement and TGDM mask routing audit — 2026-10-17

> **Atlas correction (2026-10-18):** the owner-selected AtlasEM target is now
> `ibsr11objV3`; the historical `ibsr8obj` unavailability statement below is
> superseded. Native source files are packaged, while loader adaptation,
> checked conversion, and selected-atlas validation remain open. See the
> [current atlas findings](plus-atlas-selection-2026-10-18.md).

## Final MaSeg refinement

The selected `run_maCRUISEplus_v1` route always calls `run_SegmentationRefine`
after `run_CruisePlusSurf_pipeline`. The restartable Python raw-input runner
previously stopped after surface and thickness production. It now appends a
`maseg_refinement` stage using the existing source-mapped V3 helper, with the
selected 0.2 mm coefficients, central filling, radius 2, 50-sweep cap, and
source-faithful mode. The result is saved as full-grid `SEG_refine.nii.gz`,
recorded in `final_outputs`, and returned by the workflow result. The stage
records input/output hashes and rule metrics; resume verifies the output hash
and repairs a missing or corrupted result. This closes a concrete integration
omission; it does not establish complete PLUS equivalence.

## Selected PLUS TGDM mask argument

The selected MATLAB call chain is `run_CruisePlusSurf_pipeline.m` →
`run_TGDMinvcrop.m` → `run_TGDMinvcrop_jist.m` → JIST `MaTGDMcrop`.
`run_TGDMinvcrop_jist.m` sets the initial membership, WM, GM, ACE-GM, thinned
skeleton, GVF, and crop inputs. It does not set the optional `Brain Mask Level
Set` input. In `MaTGDMcrop.java`, that optional input becomes Java `null` when
unset (lines 298–301), and that null is passed to all three selected inner,
central, and outer `MaGenericTGDM` calls. The Python full PLUS route had been
injecting a unit signed field derived from the supplied binary brain mask. It
now leaves the optional TGDM bound unset on this selected wrapper route. The
user mask remains active in the source-supported input preparation and crop
path. The lower-level Python composition API still supports an explicit bound
for callers that need that separate mode.

Source file SHA-256 values audited:

| Source | SHA-256 |
|---|---|
| `pipeline/run_CruisePlusSurf_pipeline.m` | `aadb7628e844d9795e52a184920e570a979efaf930c0e7199c034a857413e059` |
| `pipeline/run_TGDMinvcrop.m` | `b21c25b8bb4edfa8b554309a45506b4d47764c81a4b286f7c8af9b4bd5a91381` |
| `pipeline/run_TGDMinvcrop_jist.m` | `90a40f3330adab3b21615afcc5c373c44637e3a887bed9e509c83c1efc172131` |
| `toads-cruise/.../tgdm/MaTGDMcrop.java` | `1d4b5540dd5ea13c2b6eed5235c7ae05574432c39ff4240b0d595ea0342042d0` |

The 2026-10-17 [mask-aware held-fixed comparison](in-vivo-plus-tgdm-mask-aware-java-parity-2026-10-17.md)
remains valid for the explicitly supplied optional-bound mode. It is not
evidence for the selected PLUS wrapper because its Java harness deliberately
provided a bound absent from the wrapper. The earlier unmasked comparison is
also only held-fixed TGDM evidence; neither validates an independently run
complete PLUS source sequence.

## Remaining blockers

The owner-selected `ibsr11objV3` native descriptor/payloads are present, but
selected-loader descriptor adaptation, checked RAS conversion, and human
AtlasEM convergence remain open. The historical Java/MATLAB wrapper still
names `ibsr8obj`, so its resource selection must be reconciled before a
complete selected-target comparison. The original MATLAB/MIPAV/JIST runtime
and selected NiftyReg Aladin executable are unavailable here; the
ANTsPy-versus-Aladin optimizer comparison remains open. A fresh held-fixed
Java replay on the selected null-bound TGDM inputs is now recorded in the
[2026-10-18 comparison report](in-vivo-plus-tgdm-null-bound-java-parity-2026-10-18.md).
