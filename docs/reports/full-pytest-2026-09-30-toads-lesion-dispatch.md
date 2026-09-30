# Full-suite evidence: TOADS and LesionTOADS output/dispatch adapters

Date: 2026-09-30

## Full suite

- Command: `.venv/bin/python -m pytest -q`
- Result: **719 passed, 2 expected RAS reorientation warnings, 163.49 seconds**.
- Tested implementation commit: `bc8091c14df0a110b6d282759d498dd355b7a823` (`Pin folded TOADS dispatch Java reference`).
- The warnings are intentional. The Dura topology test supplies LPS topology and paint-mask volumes and verifies reorientation to RAS.
- Focused adapter/dispatch command: `.venv/bin/python -m pytest -q tests/test_toads_thalamus_output.py tests/test_toads_exact_dispatch.py tests/test_lesion_csf_distance.py` — **18 passed in 1.00 second**. P18cv and P39zzzzi were also included in the full run; their focused tests previously passed 6 and 2 cases, respectively.

## P18cv: CRUISE output maps

The selected `MultipleToadSegmentation.generateCruiseOutputs()` excerpt and Python adapter agree for CSF, GM, WM-fill, and hard-WM maps in exact and approximate modes. Committed true-RAS NIfTI fixtures cover aliases, overwrite order, exclusions, candidate slots, and hard-label behavior. Wrapper uncrop remains the separate P18cn operation. This does not validate segmentation evolution.

## P18cw: thalamus membership output

The selected `generateThalamusOutput()` excerpt and Python adapter match on exact and approximate membership output. Tests preserve exact `Thalamus` name matching, duplicate atlas-class overwrite, candidate-slot overwrite, absent-name zero fill, and wrapper uncropping through P18cn. This validates an optional output adapter, not thalamus segmentation.

## P18cx: fixed-input exact membership dispatch

The adapter routes only OBJECT/MASK and OBJECT/MASK/OUTLIER mixes to existing source-matched kernels. Both routes match pinned Java output bits and objective values; the folded asymmetric RAS fixture independently covers OBJECT/MASK. Tests check source class/update order, auxiliary OUTLIER output, and untouched six-face borders. LESION, OUTMASK, RAYLEIGH, optimized tissue types, and cross-family mixtures reject explicitly; full `computeExactMemberships()` dispatch and statistics estimation remain incomplete.

## P39zzzzi: stopping-step lesion subtype map

The selected `updateLesionClassifcation()` excerpt matches the Python output map and count, including stale-map clearing, strict thresholds, lesion-favoring subtype ties, and lesion-disabled behavior. This is one active stopping-step output operation, not full LesionTOADS.

## P39zzzzj: dormant CSF-distance helper

The selected `computeCSFDistance()` helper and existing signed-byte `fastDistance` primitive match the pinned Java excerpt on exact-name seed masks and distance maps. True-RAS cases cover duplicate and near-miss labels, empty/border/interior support, cap 1, and default cap 5. The selected wrapper comments out this call, leaves `CSFDist` zero, and does not export it; this is a dormant support method, not active LesionTOADS behavior.

Each method remains individually scoped in [plan.md](../../plan.md). Fixture provenance is in the adjacent `tests/data/synthetic` README files for each adapter.
