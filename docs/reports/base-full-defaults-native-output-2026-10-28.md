# BASE full defaults and original-grid delivery — updated October 29, 2026

The included wheel was tested in a separate installation. A complete,
independent BASE run using the full defaults finished on Adam's supplied
T1/MaSeg/brain mask. All nine stages completed; both TOADS passes reached
their energy thresholds naturally, at iterations 20 and 23. The manifest
reports `converged: true`. This describes termination, not source equivalence.

## Corrected behavior

- Full default settings: registration 4 levels/50 iterations; post-refresh/max
  registration 5/5; TOADS 50; Dura 5; ACE 25; GVF 200; initial/inner/central
  TGDM 6/4/7; WM iso-level 0.5; updated-CSF mode 2; source-faithful refinement.
- TOADS uses axial LPS processing order internally. Exported image orientation
  is recorded independently.
- N4 operates on the original acquisition before working-grid conversion and
  retains the NIfTI reader's stored spacing. The corrected adapter matches
  bare current ANTs on all 9,830,400 voxels, with zero absolute error.
- Final NIfTIs use the original T1 dimensions, affine, qform/sform and storage
  axes. Default axis-only processing is reversed without interpolation at
  final export; an explicitly regridded run is resampled back to the original.
- Labels use integer storage; final meshes use original physical RAS mm.
  Native export is a separate restartable stage and supplies `final_outputs`.

## Engineering evidence

The final focused regression selection passed 70 tests (22 expected orientation
warnings). A separate installation of the final wheel passed 27 checks (13
expected warnings). Another 13 shared N4/PLUS tests passed. These selections
are recorded separately; they are not a count of unique tests across the
entire project. The complete project test suite was not rerun for this bundle.

Tests include oblique and non-RAS geometry, integer labels, analytic physical
fields, mesh coordinates/scalars, CLI/API defaults, independent TGDM overrides,
restart/corruption recovery, N4 reader parity and native-before-working-grid
processing. The completed real scan passed the included output checker:
80 stage artifact hashes verified; all five final NIfTIs exactly match the
original 256×256×150 shape and affine; SEG_refine is int32 with zero fractional
voxels. Both TOADS convergence checks passed.

Adam's previously supplied Python SEG_refine also had integer storage and
integer values. The screenshots alone did not establish a fractional-value
bug. Viewer interpolation should be set to nearest-neighbor for label display.

## Comparison with the supplied original MaCRUISE output

No fitted registration, smoothing, label cleanup or component removal was
applied to make these results agree. Original raw VTK coordinates are voxel
indices; applying the original Raw/Seg affine agrees with its supplied FreeView
ASC coordinates to less than 0.000009 mm. The comparison uses that verified
physical frame for the original and RAS mm for Python.

The new segmentation has **97,837 differing voxels out of 9,830,400** (about
1.0% of the image), compared with 585,542 in the old short run on the original
grid. **84 of 134 nonzero labels have Dice below 0.97**. Mean per-label Dice
is 0.93822. Background dominance makes the all-image percentage insufficient
as an agreement measure.

| Surface | Symmetric mean distance (mm) | Sampled maximum distance (mm) | Python components | Original components |
|---|---:|---:|---:|---:|
| Inner | 0.36073 | 7.73516 | 36 | 1 |
| Central | 0.34652 | 8.27032 | 36 | 1 |
| Outer | 0.26585 | 8.22603 | 36 | 1 |

Distances sample every vertex against the opposite triangle mesh; the sampled
maximum is not a continuous triangle-to-triangle Hausdorff bound. All three
Python meshes are closed by edge incidence, but extra components remain.
The comparison does not pass the project gates of per-label Dice ≥0.97,
symmetric mean surface distance ≤0.01 mm and sampled maximum ≤1 mm.
The coordinate/grid correction is verified; reconstruction equivalence and
Adam's concern about surface anatomy remain open.

## Located and unresolved upstream differences

Refinement replay using the original MaSeg and original inner/outer fields
matches the supplied original SEG_refine exactly (zero differing voxels).
This isolates the selected refinement operator; it does not validate the
independent reconstruction or its upstream inputs.

Current N4 still differs from the unidentified original backend (masked
correlation 0.98743), even though the corrected Python adapter is exact against
current bare ANTs. The original container/backend version is unknown.
Source inspection also identified a windowed-sinc preprocessing step in the
original SPECTRE boundary that is not integrated into this release. A separate
prototype matched the selected Java operation on this case; it is excluded
because the complete preprocessing handoff and downstream consequences still
need validation. No claim assigns all remaining errors to either difference.

`validation/comparison.json` preserves full quantitative results and limitations.
The bundle includes no patient arrays, reconstructed patient surfaces, original
reference programs, Java/MATLAB harnesses, containers or experimental prototype.
`SHA256SUMS` covers every shipped file except itself. Upstream license/notice
texts are preserved; source code and CVS metadata in the notice tree are excluded.

## Local evidence and deliverable

- Completed run: `build/adam-full-native-n4-2026-10-28/run.json`
- Detailed comparison: `build/adam-output-triage-2026-10-28/full-defaults-comparison.json`
- Visual comparison: `build/adam-output-triage-2026-10-28/full-profile-native-grid-comparison.png`
- Python-only handoff: `build/pymacruise-adam-0.1.1-2026-10-29.tar.gz`
- Source implementation signature: `4e47d088aba9ab612fbadc56db862d0b196235d15f2ef8e223186c3201728ed2`

The received folder remains unchanged. The earlier diagnostic runs were interrupted
when adapter/order differences were found; only the native-N4 run above is the
completed final run used in these comparisons. Bundle-relative validation links
in the preceding sections refer to the extracted handoff directory.
