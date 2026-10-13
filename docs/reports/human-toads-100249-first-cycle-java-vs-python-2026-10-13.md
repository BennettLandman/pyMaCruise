# 100249 human TOADS first-cycle Java/Python comparison (2026-10-13)

## Scope

This is a held-fixed, single-cycle comparison of the selected Java VandyCRUISE TOADS wrapper and Python on case 100249's native RAS T1/mask grid. It isolates the first post-registration classification update; it does not establish full-loop or end-to-end MaCRUISE equivalence.

The input is the saved SPECTRE-preprocessed `stripped_t1-RAS.nii.gz` (176×240×256; spacing 0.9999984, 1, 1 mm) and matching all-voxel support mask. Java receives the source x-fastest float32 serialization in `build/human-toads-source-pipeline-100249/java-input.raw`; Python reads the same RAS images. The atlas is the bundled 10-object VandyCRUISE descriptor. Java and Python output use the same 11 labels and grid.

## Method

- Java ran the source harness with `maximum_iterations=1`, allowing its first growth and final label update. The raw Java result is SHA-256 `311a01d53acb18ee20a59c66aeec2ca560e29c05f807cad4165185fddd43a7ed`.
- Python ran the source-selected registered startup, exact first membership pass, first registered loop registration, first growth, and source label update. The Python run was capped at the next iteration boundary; it did not continue the full TOADS loop. Peak sampled RSS was 5.120 GB under an 8 GB watchdog; no RSS termination occurred.
- Java's MIPAV x-fastest map was reshaped in Fortran order and compared voxel-by-voxel with the Python RAS NIfTI.

## Results

There are 5,445 differing voxels of 10,813,440 (0.05035%). Every label meets the agreed minimum Dice of 0.97; the minimum is label 17 at 0.98110. Label 22 is 0.98324 and label 6 is 0.98466. All 11 labels pass this threshold at this single-cycle checkpoint.

Of the 5,445 differences, 5,173 lie on a 26-neighbor boundary in the Python map and 5,420 lie within one voxel of that boundary. This supports a boundary-localized mismatch; it does not identify the cause. Between cycle-one and cycle-two outputs, 3,863 cycle-one disagreements resolve, 1,582 persist, and 16,975 new disagreements appear; Java changes 324,708 labels and Python changes 323,787 labels between those checkpoints. These transition counts are in [cycle-transition metrics](human-toads-100249-cycle-transition-metrics-2026-10-13.json). The Python NIfTI SHA-256 is `618d74a2dc324841a6885fa0b77adbb5385ed9ba045f3afab91f45e1afcf82b4`.

Detailed per-label counts and Dice values are in [per-label JSON comparison metrics](human-toads-100249-first-cycle-java-vs-python-2026-10-13.json). Boundary-localization counts are in [boundary metrics](human-toads-100249-first-cycle-boundary-metrics-2026-10-13.json). The full images and raw Java buffer remain in ignored `build/`; input identity is retained in the preceding reports and the local `java-input.json` run record.

## Interpretation and remaining work

This first-cycle map meets the user's 3% segmentation tolerance for each TOADS label on this held-fixed human input. It substantially narrows the prior two-iteration mismatch: that comparison had 18,557 differences and label 17 Dice 0.9556. The first-cycle disagreement alone does not explain the cycle-two failure: the maps evolve differently enough that 16,975 new disagreement voxels appear by cycle two. The later-cycle divergence therefore remains to be localized; the natural-stop and historical LPS/RAS discrepancies also remain open. Do not mark P18cycaq validated based on this one cycle.

The result is source-comparison evidence for a selected TOADS stage on one human input, not an independent anatomy validation, not proof of all TOADS behavior, and not whole-workflow equivalence. The native RAS output contract remains in force.
