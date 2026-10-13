# Superseded atlas-aligned synthetic PLUS baseline — 2026-10-12

**This initial report used an incorrect synthetic MaSeg label mapping and is retained only as a historical record. Its ribbon metrics are not evidence for the corrected fixture.** The old generator grouped topology classes 14, 16–18, and 22 as cortical GM and class 24 as cortical WM. The source-derived mapping is now corrected and regression-tested: class 15 alone maps to GM 100, classes 25–27 map to WM 44, classes 5–6 map to CSF 5, class 1 maps to background 0, and remaining listed structures map to other tissue 75.

The corrected generator rewrote `maseg.nii.gz` and `provenance.json`; T1, FLAIR, and brain mask remain unchanged. The original all-nine-stage run, output hash audit, and ribbon numbers documented before this correction describe the old labels only.

Use [the corrected terminal synthetic PLUS report](plus-atlas-consistent-corrected-terminal-2026-10-12.md) for current results. It records the complete corrected workflow, 95/95 manifest hash verification, 69/69 finite RAS NIfTIs, natural convergence at iteration 10, the explicit ribbon-quality failures, and the selected JIST surface extraction comparison on the identical final level sets.
