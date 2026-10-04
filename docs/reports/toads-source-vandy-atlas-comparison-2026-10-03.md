# TOADS source-atlas comparison — 2026-10-03

## Scope

This is a synthetic, one-iteration comparison between the supplied VandyCRUISE
Java implementation and PyMaCRUISE. Java used the original 11-structure raw
atlas descriptor and Python used the committed converted NIfTI atlas. Both ran
against the same `(36,45,34)` synthetic T1 values. The Java input was exported
in Fortran order to preserve MIPAV's x-fastest voxel convention. This is not a
whole-workflow equivalence result and is not real-data validation.

## Atlas conversion check

The conversion audit and source-payload round-trip check covered 35 atlas
images: 32 raw payloads restored exactly using their source dimensions,
datatype, and MIPAV voxel order; the three Analyze images were voxel-exact when
canonicalized through their source header affine. The Analyze source-axis
metadata still disagrees with the conversion audit (`L,A,S` in the headers vs
`L,P,S` in the audit), so physical origin/direction claims for those examples
remain unverified. Those three SPECTRE images are not used in this Vandy atlas
comparison.

## Java/Python comparison

Java completed one source outer iteration using `VandyCoreSyntheticRunExport`
with a 3 GB maximum Java heap. Python completed the registered trace with one
post-refresh registration update and two iterations (the Python loop requires
at least two). Both runs used one registration level, one registration update
per branch, and relation limit 3.

Java's startup membership and transformed-prior grids are `(53,60,51,11)`;
Python's are `(52,60,51,11)`. The Java debug log reports crop bounds
`[3,35] [3,42] [0,30]`; Python reports `[2,33] [3,42] [0,30]`. With ten
padding voxels, their shared x extent maps Java planes `0:51` to Python planes
`1:52`. On that overlapping region:

- Startup memberships differ by mean absolute `0.0073713`, maximum `1.0`, and
  124,459 of 1,716,660 values above `1e-6` (7.25%).
- Transformed shape priors differ by mean absolute `0.0098571`, maximum
  `1.0`, and 74,023 of 1,716,660 values above `1e-6` (4.31%).
- The categorical transformed topology agrees exactly at 143,642 of 156,060
  voxels (92.04%).

These are comparisons on the crop-derived overlapping grid, not an equivalence
result. The source and Python startup objective values also differ, but are
not treated as comparable until the crop and initialization are aligned.

The Python run reports iteration energy/ratio `105.751678 / 0.70089215` at
iteration 1; Java reports `50.30009 / 1.0909222`. These values confirm a
substantial downstream difference but do not localize it while startup grids
differ. No attempt is made to tune either implementation to these totals.

## Evidence and next action

- Java source log: ignored diagnostic file `build/java-reference/source-vandy-atlas-first-loop.log`; source completed with output dimensions `[36,45,34]`.
- Python source-atlas log: ignored diagnostic file `build/java-reference/python-source-vandy-first-loop.log`; run completed in 19.7 seconds and produced RAS uint8 labels.
- Synthetic inputs: committed fixture `tests/data/synthetic/toads_downsampled_vandy_atlas/` and bundled atlas `data/atlases/cruise-atlas-10obj-toads2008.txt`.
- The captured arrays are ignored build artifacts, not committed as test samples.
- No pytest suite was run for this documentation-only evidence update.

Next, explain the one-voxel startup crop difference and compare constructor
labels, centroids, variances, and memberships after the bounds and coordinate
mapping match. Do not infer full TOADS equivalence or workflow readiness from
this partial overlap comparison.
