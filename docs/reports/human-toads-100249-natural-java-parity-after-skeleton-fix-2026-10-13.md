# 100249 natural-convergence TOADS Java/Python parity after skeleton fix (2026-10-13)

## Result

The corrected Python TOADS loop matches the selected Java RAS natural run exactly on the held-fixed case 100249 T1, all-voxel mask, VandyCRUISE descriptor, and source settings. The terminal classifications have **0 differing voxels of 10,813,440**; Dice is **1.0 for all 11 labels**. Both naturally stop after 38 iterations. Python's terminal energy ratio is `1.3296979e-05`, below the `1e-4` threshold. The Python output is RAS.

This follows the fix that forwards the Java `MultipleToadSegmentation` constructor skeleton into first-cycle thinning. The constructor seeds border voxels from the initial classification and the first label initialization rewrites only interior voxels. Corrected matched-input comparisons also match exactly at one and two cycles; see the [startup/skeleton investigation](human-toads-100249-cycle-replay-after-skeleton-fix-2026-10-13.md).

## Reproduction and resource use

The natural run used `pymacruise toads-loop` with 4 registration levels, 50 registration iterations, 5 post-refresh iterations, 5 registration maximum iterations, and maximum distance `1e-4`; no iteration cap override was supplied. The input grid is `176×240×256` (10,813,440 voxels). Java raw data was compared in its x-fastest Fortran reshape order against the Python RAS NIfTI.

Peak sampled Python RSS was **6.488 GB** under an **8 GB watchdog**; the run exited normally without watchdog termination. Detailed label counts and exact voxel agreement are in the [machine-readable comparison](human-toads-100249-natural-java-parity-after-skeleton-fix-2026-10-13.json).

## Scope and remaining gap

This is source-comparison evidence for the natural RAS TOADS loop on one held-fixed human input. It does not show that the historical LPS preprocessing path can be converted to the same RAS result: the saved first-pass labels match Java in the historical LPS orientation, while Java's LPS and RAS runs already differ during crop/registration. That orientation-sensitive source behavior remains open for P18cycaq. This result is not independent anatomical validation or whole-workflow MaCRUISE equivalence.
