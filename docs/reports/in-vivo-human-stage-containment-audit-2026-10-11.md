# Human MaSeg Stage-Containment Audit

> **Historical, preliminary audit.** Its source run selected ACE mode 2 but
> omitted the V3 cleaned MaSeg boundary input. Use the corrected metrics and
> boundary evidence in [the boundary-enabled human tail report](in-vivo-base-tail-with-v3-boundary-2026-10-11.md)
> and its ignored `build/invivo-base-source-v3-human-boundary-2026-10-11/`
> outputs for current human tail results.

Evidence: **read-only in-vivo stage-boundary localization; not source equivalence, clinical validation, or causal proof**.

Run: `build/invivo-base-source-v3-human-rasgrid-2026-10-11`; status `complete`.

Raw input MaSeg counts: GM 642,291; WM 428,303.

Support uses raw input labels, not a source reference. Continuous maps report fraction with value >102 (0.4×255); binary masks use >0. TOADS categorical class images are inventoried but excluded because class indices were not mapped to MaSeg GM/WM labels. TGDM GM is inner-positive/outer-nonpositive ribbon and WM is inside inner (≤0). These are containment diagnostics, not Dice scores.

| Saved checkpoint | GM support | WM support | Domain / interpretation |
|---|---:|---:|---|
| TOADS first GM membership | 90.65% (582,244/642,291) | 15.81% (67,694/428,303) | full input grid; scaled membership > 0.4*255 |
| TOADS first WM-fill membership | 4.61% (29,624/642,291) | 89.06% (381,433/428,303) | full input grid; scaled membership > 0.4*255 |
| TOADS first WM mask membership | 2.69% (17,285/642,291) | 83.68% (358,412/428,303) | full input grid; nonzero membership |
| TOADS second GM membership | 90.98% (584,332/642,291) | 15.68% (67,167/428,303) | full input grid; scaled membership > 0.4*255 |
| TOADS second WM-fill membership | 4.87% (31,305/642,291) | 89.19% (382,015/428,303) | full input grid; scaled membership > 0.4*255 |
| TOADS second WM mask membership | 2.84% (18,215/642,291) | 83.83% (359,064/428,303) | full input grid; nonzero membership |
| Dura corrected GM map | 88.79% (570,291/642,291) | 15.79% (67,642/428,303) | full input grid; scaled membership > 0.4*255 |
| Dura corrected WM map | 4.61% (29,604/642,291) | 89.05% (381,420/428,303) | full input grid; scaled membership > 0.4*255 |
| Dura GM topology output | 97.32% (625,058/642,291) | 99.99% (428,260/428,303) | full input grid; nonzero; not treated as probability |
| Dura WM topology output | 42.59% (273,560/642,291) | 99.91% (427,922/428,303) | full input grid; nonzero; not treated as probability |
| Dura brain mask | 96.76% (621,452/642,291) | 99.98% (428,236/428,303) | full input grid; binary nonzero mask |
| TOADS GM transferred to input grid | 90.98% (584,332/642,291) | 15.68% (67,167/428,303) | full input grid; scaled membership > 0.4*255 |
| TOADS WM-fill transferred to input grid | 4.87% (31,305/642,291) | 89.19% (382,015/428,303) | full input grid; scaled membership > 0.4*255 |
| Fused GM in Cruise crop | 92.03% (591,079/642,291) | 12.90% (55,242/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| Fused WM-fill in Cruise crop | 4.87% (31,305/642,291) | 90.76% (388,742/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| Fused mask in Cruise crop | 2.84% (18,215/642,291) | 86.93% (372,312/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; nonzero mask |
| TOADS GM in Cruise crop | 90.98% (584,332/642,291) | 15.68% (67,167/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| TOADS WM-fill in Cruise crop | 4.87% (31,305/642,291) | 89.19% (382,015/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| Cruise input GM crop | 92.03% (591,079/642,291) | 12.90% (55,242/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| Cruise input WM-fill crop | 4.87% (31,305/642,291) | 90.76% (388,742/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| Cruise input WM mask crop | 2.84% (18,215/642,291) | 86.93% (372,312/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; nonzero mask |
| Cruise topology-corrected WM fill | 5.22% (33,541/642,291) | 90.72% (388,542/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| ACE enhanced GM | 80.56% (517,416/642,291) | 12.89% (55,188/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; scaled membership > 0.4*255 |
| ACE skeleton | 12.78% (82,095/642,291) | 0.01% (61/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; nonzero occupancy |
| ACE thinned skeleton | 7.80% (50,077/642,291) | 0.01% (38/428,303) | voxel-aligned crop of input grid; crop start [30, 58, 105]; nonzero occupancy |
| TGDM terminal GM containment | 82.34% (528,870/642,291) | — | full input grid; GM labels in inner>0 and outer<=0 ribbon |
| TGDM terminal WM containment | — | 88.48% (378,979/428,303) | full input grid; WM labels inside inner<=0 |
| Refined segmentation GM containment (denominator-only comparison) | 88.75% (528,870/595,895) | — | full input grid; same terminal TGDM region evaluated on refined labels |
| Refined segmentation WM containment (denominator-only comparison) | — | 96.48% (378,979/392,800) | full input grid; same terminal TGDM region evaluated on refined labels |

## What the saved boundaries show

Across TOADS GM and WM-fill maps, input-label membership values remain tissue-selective at the >102 screening threshold. After fusion into the Cruise crop, GM-label overlap with high GM membership is 92.03%; ACE-enhanced GM is 80.56%. WM-label overlap with high WM-fill is 90.72% in the fused crop and 90.72% after WM topology correction. Terminal TGDM containment is 82.34% for GM ribbon and 88.48% for WM. This localizes the first clear saved-boundary decrease in high-GM membership to ACE enhancement, while terminal containment is measured after TGDM. These are different metrics (map threshold vs geometric containment), so this pattern is a clue, not causal proof.

Dura topology outputs are recorded as nonzero occupancy only. Their numeric values are not interpreted as probabilities. The Dura brain mask contains 96.76% of raw GM labels and 99.98% of WM labels; this differs from supplied-mask vs Dura-mask Dice and does not establish source equivalence.

## Localization limits

The run saves TOADS first/second classifications and tissue maps, Dura maps, transferred/fused maps, Cruise crop and derived ACE maps, and terminal TGDM fields. It does **not** save per-iteration or per-sub-evolution TGDM level sets, so exact loss point inside TGDM cannot be identified. Terminal fields show support by the end of Cruise core; intermediate membership maps are screening trends and do not prove causal attribution.

TGDM iteration checkpoint search found 0 matching NIfTIs. The manifest’s legacy `tgdm_iterations=1` is superseded by explicit initial/inner/central options 6/4/7; core provenance records the wrapper-mapped outer count as 7. Per-step TGDM sub-evolution fields are not saved. See JSON for all metrics, paths, label IDs, and crop geometry.

Machine-readable metrics, paths, label IDs, and crop geometry: `build/invivo-base-source-v3-human-rasgrid-2026-10-11/stage-containment-audit/stage-containment-audit.json` (ignored local build output).

This audit does not demonstrate source equivalence, clinical validity, or causal proof.
