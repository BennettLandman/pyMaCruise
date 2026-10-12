# P41z iteration-15 downstream ribbon audit — 2026-10-11

## Run and result

The iteration-15 nonterminal registered-loop checkpoint was passed through
the source-scaled LesionTOADS Dura wrapper, native transfer, and plus
topology/ACE/GVF/TGDM/surface tail. The run completed under the 4 GB RSS
watchdog at **0.672 GB peak**. The resulting surface manifest is
`build/p41z-case-normalization/iterations-0012-to-0015-surfaces/plus-surfaces/run.json`.

| Measure | Iteration 11, corrected wrapper | Iteration 15 | Acceptance |
| --- | ---: | ---: | ---: |
| WM label 44 inside inner | 99.406% | 99.435% | ≥97% |
| GM label 100 between inner/outer | 92.835% | 92.890% | ≥97% |
| CSF label 5 outside outer | 99.258% | 99.258% | ≥97% |
| Closed connected mesh components | 2 per surface | 1 per surface | 1 |
| Mesh Euler characteristic | −12 per surface | −18 per surface | 2 |
| Surface-thickness median | 6.434 mm | 6.496 mm | 4–6 mm |
| Surface-thickness p05 / p95 | 0.237 / 12.280 mm | 0.236 / 12.298 mm | 2–8 mm |
| Level-set mean thickness | 4.832 mm | 4.842 mm | 4–6 mm |
| Nested voxel regions | Pass | Pass | Pass |

The extra four registered-loop passes slightly improve GM and WM support and
remove the extra connected component, but they do not bring GM support or
surface-thickness distribution within bounds. The surfaces remain topologically
wrong despite being closed and connected; Euler −18 does not meet the required
genus-zero value of 2. Mean level-set thickness passes, but this does not erase
the surface median and percentile failures.

## Evidence scope

This is a held-fixed synthetic P41z continuation. Iteration 15 is still
nonterminal at energy ratio **0.13346674**, far above the selected `1e-4`
threshold. The continuation ratios through iterations 8–15 and timing are in
the [registered-loop report](p41z-registered-loop-iterations-12-to-15-2026-10-11.md).
The downstream Java TGDM comparison passes its selected surface gates, as
recorded in the [corrected-input TGDM report](p41z-scaled-wrapper-tgdm-java-2026-10-11.md),
but the current run does not prove whole-loop equivalence or real-human-MRI
acceptance.
