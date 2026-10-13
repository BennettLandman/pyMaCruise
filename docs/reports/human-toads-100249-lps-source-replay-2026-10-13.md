# 100249 saved TOADS source replay and RAS mismatch — 2026-10-13

## Finding

Case 100249 provides a stronger saved reference for the historical TOADS input
than the earlier human sample. The Java source harness reproduces the saved
first-pass classification **exactly** when it reads the saved
`T1_N4Corrected_MAstripped.nii.gz` array in its original LPS axis order and is
given the corresponding MIPAV axial LPS axis directions (`R2L`, `A2P`, `I2S`).
After a simple x/y axis flip to RAS array order, Java and the saved
`Toads1stPass/*toads1seg.nii.gz` have 0 changed voxels of 10,813,440 and Dice
1.0 for all 11 labels.

The Python RAS run does **not** match that saved reference. Python naturally
stopped at iteration 33 with terminal energy ratio `5.25178002e-05`; Java's
LPS replay also produced 33 iteration records, but its map differs from
Python at **239,301 voxels** after the x/y orientation change. Per-label Dice
for Python versus Java LPS ranges from **0.101 to 0.998**; labels 16, 17, 22,
and 24 are below 0.65. Thus the iteration count is not the cause by itself.

A Java replay that instead labels the already-RAS array as RAS differs from
Python at 66,936 voxels (per-label Dice 0.879–0.999). A 33-iteration-capped
Java RAS replay differs from Python at 67,135 voxels. These runs establish
that the current selected Python implementation does not reproduce Java on
this source-derived intensity under the RAS path, even at a matched iteration
cap. This contrasts with the exact two-pass comparison on the earlier human
case in [the matched TOADS report](human-toads-matched-java-comparison-2026-10-13.md).

## Inputs and comparison boundaries

- Historical TOADS input: `100249/outputs/ToadsSeperate/MultiAtlasStripping/T1_N4Corrected_MAstripped.nii.gz`; shape `176×240×256`, 1 mm, LPS; SHA-256 `7eb0ae628b618454cb26211d35f12db70e6ea1bbf3cee57dd9a30b2f770bd126`.
- Saved source first-pass labels: `100249/outputs/ToadsSeperate/Toads1stPass/T1_N4Corrected_MAstripped_toads1seg.nii.gz`; shape `176×240×256`, LPS orientation with zero-origin header; SHA-256 `bdc3390191a01379f02a757c15eb1953b5b01cbca2f784e11d2fb33cc802450e`.
- The Java LPS replay used the same float32 image voxels, an all-voxel support mask equivalent to the already-stripped source image, the bundled selected descriptor, and the selected VandyCRUISE settings (4 levels, 50 registration iterations, 5 post-refresh iterations, 50 maximum segmentation iterations, `1e-4` convergence and maximum-distance, relation limit 3). It completed after 33 recorded iterations; the harness does not expose a terminal stop reason.
- The Python run read the physically matching image after reorientation to RAS and used an all-voxel support mask. Its public output is RAS. Comparison to the saved source labels is by x/y-flipped voxel indices only: saved TOADS headers have a zero origin and do not establish physical-world alignment with the SPECTRE-preprocessed T1.
- The Java LPS output after x/y reorientation has SHA-256 `efb425a117df2b9aaa7d974ec06161c9c099550f17e8ad887340ed662ed1f3af`; the Python RAS output has SHA-256 `410972a9868ba2d4d2db173ce88ec19a223ac160d30d9d3860e22eea30266f1c`.
- Peak RSS was 3.744 GB for the Java LPS replay and 7.466 GB for Python, under separate 16 GB watchdogs.

The Java test harness now accepts an optional `RAS` or `LPS` orientation
argument while retaining `RAS` as its default. It was compiled and run against
the VandyCRUISE source checkout and the selected atlas descriptor. The source
image, saved source outputs, and all replay intermediates remain outside the
repository or in ignored `build/`; none are committed.

## Interpretation and next work

This run demonstrates that the historical saved first-pass map can be
reproduced by the selected Java implementation when the original LPS array
and axis directions are used. It also shows that reorienting that image to
RAS before the current Python TOADS path changes the result materially. The
source method therefore has an orientation-sensitive behavior that the
current canonical-RAS path does not yet preserve on this case. The result is
an algorithm/source integration failure to resolve, not a validation pass.

The checked source locates one concrete place where this information matters:
`MultipleToadPreprocess` sets axial axis signs from the image's original
`FileInfoBase` directions, and `getSignedImageResolutions()` divides voxel
spacing by those signs. Thus the selected source supplies
`(+sx,+sy,+sz)` for an LPS-ordered array and `(-sx,-sy,+sz)` for a RAS-ordered
array. Python currently normalizes the voxel data to RAS first and then always
uses the latter signs. This is a verified source-metadata difference, but it
does not yet prove that it alone causes the output mismatch. The next useful
comparison is to expose matched initialization, crop, and registration states
before changing the workflow. A potential compatibility adapter can retain
the original array order internally and convert segmentation outputs to RAS;
that design still needs source-state evidence and must keep every saved output
in RAS.

The Java logs provide an early localization lead: its LPS and RAS replays
produce different post-crop translations (`[-0.4999992, 6.5, 14.0]` versus
`[-0.4999992, 6.0, 14.0]`) and different converged registration transforms.
This shows divergence is present by registration, before the label-evolution
loop. The current Python human run did not retain equivalent startup telemetry,
so the next comparison should record Python crop bounds, transformed atlas,
initial memberships, and registration transform against both Java runs.

Next, compare the Java LPS and RAS registration/crop/initial-membership states
with the Python equivalents and determine where the axis reorientation
changes the branch. Any correction must preserve the project contract that
saved outputs are RAS, match all labels to Dice ≥0.97, and keep the input
orientation warning/provenance. Then resume only downstream BASE stages from
the validated checkpoints. The earlier exact first/second-pass comparison
remains valid for its own held-fixed input pair; it does not generalize to
this source-preprocessed LPS case.

Machine-readable comparisons:

- [Java LPS versus saved source labels](human-toads-100249-lps-vs-source-2026-10-13.json)
- [Java LPS versus Python RAS](human-toads-100249-lps-vs-python-2026-10-13.json)
- [Java RAS versus Python RAS](human-toads-100249-ras-java-vs-python-2026-10-13.json)
