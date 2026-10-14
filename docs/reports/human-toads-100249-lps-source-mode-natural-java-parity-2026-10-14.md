# 100249 source-axis-compatible TOADS natural parity (2026-10-14)

## Result

An explicit source-axis compatibility option now reproduces the historical LPS VandyCRUISE run while writing Python outputs in RAS. On case 100249, the Python natural run stops at iteration 33 with energy ratio `4.63042579e-06` (threshold `1e-4`). Its terminal RAS labels match both the selected Java natural LPS output reoriented to RAS and the saved historical `Toads1stPass` labels reoriented to RAS exactly: **0 changed voxels of 10,813,440; Dice 1.0 for all 11 labels**.

The option `pymacruise toads-loop --preserve-input-axis-order` first validates and warns about non-RAS inputs, reconstructs the source's axial RAS/LPS voxel order for TOADS internals, then reorients the classification and flattened Dura wrapper maps back to RAS before writing. The default path remains canonical RAS. Synthetic CLI validation confirms the classification plus all four Dura products are RAS NIfTIs.

## Method and resources

The source T1 and same-grid all-voxel brain mask were LPS (`L`,`P`,`S`). Java used `R2L`, `A2P`, `I2S` axis metadata, with 4 registration levels, 50 registration iterations, 5 post-refresh iterations, maximum distance and convergence threshold `1e-4`, and the bundled 10-object VandyCRUISE descriptor. Python used the same settings and natural stop behavior. The human run used the 50-iteration maximum; it stopped naturally at 33. Python output orientation and geometry were checked as RAS.

Python peak sampled RSS was **7.739 GB** under an **8 GB watchdog**, and it exited normally. Focused validation: 40 image/geometry/alignment/output tests passed in 2.13 seconds; the source-axis end-to-end CLI test passed in 1.61 seconds. The full-suite evidence predates this implementation and is not claimed for it here.

Machine-readable Java/Python and saved-source comparisons are in the [comparison metrics](human-toads-100249-lps-source-mode-natural-java-parity-2026-10-14.json). Raw images and Java buffers remain in ignored `build/`; no subject images are committed.

## Scope

This closes the 100249 selected TOADS orientation compatibility comparison for this LPS case. The compatibility option currently accepts axial RAS and LPS input orders; other permutations require a separate source-validated mapping. It does not validate independent anatomy, the complete BASE/PLUS workflow, or other acquisition preprocessing chains.
