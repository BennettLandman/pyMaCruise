# 100249 two-iteration RAS TOADS comparison — 2026-10-13

## Result

The selected Java wrapper and Python implementation were run from the same
RAS T1 array, same all-voxel support mask, bundled VandyCRUISE descriptor, and
same wrapper settings. Both were capped after two segmentation iterations
(so this is a held-fixed early-iteration comparison, not a natural-convergence
run). The terminal maps differ at **18,557 of 10,813,440 voxels** (0.172%).
Label Dice values range from 0.9556 to 0.9997; label 17 is below the agreed
0.97 criterion. All other labels meet 0.97 on this two-iteration check.

The RAS Python startup probe in the [startup record](human-toads-100249-python-startup-2026-10-13.json)
matches Java's logged initial crop transform and first registration transform
to printed precision. The difference therefore appears after initial atlas
registration, within post-refresh initialization or the registered loop.
The full natural-run RAS comparison differs at 66,936 voxels, so this early
comparison establishes that the terminal mismatch is not caused only by a
different natural stopping iteration.

## Reproduction and limits

- Input grid: `176×240×256`, approximately 1 mm, RAS; masked T1 SHA-256
  `26412b0939d3303ac088330f772196f08bf9461342f94e7147a937ee0ad69291`.
- All-voxel mask SHA-256:
  `79c7775fcb4132303dd1d32e8c3b54273c6dad98c2a7383a5f5e084b02cc48bf`.
- Java two-iteration output SHA-256:
  `483a1abba571ce5dc56fcd6fb1441a71b10c498250144a4752950d5153f15a22`.
- Python output was produced by `pymacruise toads-loop --maximum-iterations 2`.
- Peak sampled RSS: Java **3.956 GB**, Python **5.974 GB**, each under an
  8 GB watchdog.
- Per-label voxel counts and Dice values are in the [machine-readable
  comparison](human-toads-100249-two-iteration-java-vs-python-2026-10-13.json).

This narrows the 100249 investigation but does not validate TOADS equivalence.
The Java LPS replay still reproduces the saved first-pass map exactly, while
LPS and RAS source runs differ at crop/registration. The RAS two-iteration
run is also below threshold for label 17. The next diagnostic is to compare
the first registered cycle's initialization and outputs, then identify the
first source-order/state difference.
