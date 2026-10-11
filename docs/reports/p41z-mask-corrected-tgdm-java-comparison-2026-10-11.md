# Corrected-mask P41z TGDM comparison — 2026-10-11

This is a selected Java/Python TGDM comparison on the corrected P41z iteration-6
crop (73×91×63, RAS), with the brain-tissue fusion mask corrected upstream.
It is a held-fixed method comparison, not a complete independent
MaCRUISEplus workflow comparison. The Java run used the selected
`MaGenericTGDM` sequence and completed under an 8 GB RSS watchdog (0.615 GB
sampled peak). The machine-readable measurements are in
[`p41z-mask-corrected-tgdm-java-comparison-2026-10-11.json`](p41z-mask-corrected-tgdm-java-comparison-2026-10-11.json).

The initial level-set output has no sign changes (maximum absolute voxel error
`4.77e-7`). The smoothed initial surface passes the surface gates (Dice 1.0,
Hausdorff `8.35e-5 mm`, symmetric mean distance `1.33e-7 mm`). Every later
surface (`initial_inner`, `inner`, `central`, and `outer`) has exactly two
sign-changing voxels and fails the 1 mm Hausdorff gate: the measured values are
`32.39–32.40 mm`. Their Dice scores are about `0.999956`; symmetric mean surface
distances are `0.006996–0.007316 mm`, within the 0.01 mm gate. Those passing
metrics do not compensate for the Hausdorff failures.

The two sign changes occur at crop indices `(26, 27, 1)` and `(34, 35, 0)`,
RAS positions approximately `(-200.391, 67.891, -240.489) mm` and
`(-184.359, 84.016, -242.499) mm`. At both positions Python is slightly
negative (`-0.0038` and `-0.0345` on the initial inner field), while Java is
positive (`16.66` and `12.34`). The first mismatch therefore appears during
inner-surface evolution, not seed generation or initial smoothing. The TGDM
cause remains unresolved; crop-edge behavior is being investigated, and no
algorithm change is justified by this comparison alone.

All six saved Python TGDM outputs were replayed exactly from the corrected
workflow inputs before comparison. The workflow remains nonterminal and its
composed ribbon-quality checks still fail. Earlier TGDM passes on the
pre-correction WM-only fusion input are not evidence for this corrected crop.
