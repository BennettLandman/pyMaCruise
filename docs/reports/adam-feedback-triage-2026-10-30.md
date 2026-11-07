# Adam's 0.1.6 feedback — October 30, 2026

The feedback confirms successful execution and native-grid export, with unresolved
anatomical disagreement and a first-pass TOADS convergence failure. A geometry
review and diagnostic tools kit is available for the existing run; reconstruction
does not need to be repeated to use it. The scientific pipeline remains 0.1.6.

## Evidence in the returned archive

The archive contains 21 JSON/log files, with no images, meshes or screenshots.
The pasted email contains image placeholders; its new visual examples were not
available for inspection. Its package/source signature exactly matches the
delivered 0.1.6 release. T1, MaSeg, brain mask and atlas descriptor file hashes
match the local test case. Scientific dependency versions also match; Adam used
Linux x86-64 and Python 3.13.14, versus the local macOS ARM64/Python 3.13.17 run.

All 84 recorded stage artifacts passed integrity checks. Final volumes have
shape 256 × 256 × 150 on the original T1 grid; SEG_refine is int32. The affine
recorded in Adam's native export equals the supplied original T1 affine exactly.

The first TOADS pass reached its 50-iteration cap without the energy threshold:
terminal ratio 0.00175488892, threshold 0.0001. The second pass converged at
iteration 16. Consequently `run.json` has `status: complete` and
`converged: false`; the verification result is a failure. The local first pass
converged at iteration 19. The Linux and macOS energy histories already differ
at the first registered iteration (0.10446737 versus 0.0989699438).

## Affine and surface interpretation

For the matching original reference supplied locally, original SEG_refine uses
MaSeg's header, while Python exports use the original T1 header, as requested.
Their maximum affine coefficient difference is 5.96e-8; the maximum physical
displacement over image corners is **0.0000152 mm**. qform codes differ (2 versus
1). This rounding-level displacement cannot explain a noticeable superior-brain
discrepancy. The original SegRefine header itself is not included in Adam's ZIP;
the new review report will record both headers directly on his machine.

Original raw MaCRUISE VTK vertices are voxel indices; Python final VTKs are
physical RAS millimetres. The comparison tool mapped original vertices through
Raw/Seg's affine and independently verified the mapping against original ASC
exports (maximum error approximately 8.6e-6 mm). A raw-VTK side-by-side view
therefore needs coordinate conversion. The exact viewer and filenames in the
new screenshots are unknown. Coordinate conversion does not eliminate the
measured anatomical differences below.

| Returned measurement | Result |
|---|---:|
| Differing refined-label voxels | 121,848 / 9,830,400 |
| Mean nonzero-label Dice | 0.92317 |
| Labels with Dice below 0.97 | 72 |
| Inner symmetric mean vertex-to-triangle distance | 0.55391 mm |
| Central symmetric mean vertex-to-triangle distance | 0.43357 mm |
| Outer symmetric mean vertex-to-triangle distance | 0.19664 mm |
| Central sampled maximum distance | 10.62232 mm |

Each Python mesh has one closed component, zero open edges and zero nonmanifold
edges. Project label and distance agreement limits fail. The sampled maximum
is not a continuous triangle-to-triangle Hausdorff bound. Good topology or
execution does not establish scientific or clinical equivalence.

## Changes made after inspecting the feedback

- Added `make-review-outputs.py`: produces explicitly paired original/Python
  labels with the same T1 header and paired meshes in physical RAS millimetres.
  Label values and mesh faces are preserved. It checks original input identity,
  rejects real grid shifts, and validates the original VTK-to-ASC convention.
  It never rewrites reconstruction products or fits a registration.
- Added both affine matrices, spatial form codes and maximum corner displacement
  to the numerical comparison report.
- Added aggregate NIfTI statistics and canonical decompressed voxel hashes to
  feedback. They distinguish image-value differences from gzip/header changes;
  images, voxel samples and mesh arrays remain excluded from feedback.
- Distinguished pipeline completion, verification success and overall exit code
  in the launcher summary. The previous summary's `returncode: 0` described the
  pipeline child; verification returned 1 and the launcher already returned a
  failure. The revised summary records both results explicitly.
- Used explicit affine arithmetic in the comparison tool and rejected nonfinite
  frame-validation errors.

Fourteen targeted helper tests pass. The paired review utility ran on the
completed local full case; saved Python mesh vertices and faces are unchanged.
Both label arrays are checked after writing. The refreshed feedback collector
was also exercised on the completed local case. These helper checks do not
validate a numerical reconstruction change; none was made.

## Remaining numerical investigation

The first differing stored scientific artifact is N4, followed by masked axial
preprocessing. The all-ones TOADS mask file is identical. Compressed NIfTI hashes
alone cannot establish the magnitude or nature of voxel differences. Matching
package versions also do not establish equality of platform-specific native
binaries. The archive lacks the corrected N4 image, so N4/backend arithmetic
versus later compiled numerical differences have not been isolated.

The concrete next step is to use the review kit on Adam's existing outputs,
return its refreshed feedback archive, and obtain
`stages/n4/attempt-01/t1-n4.nii.gz` from that same run for a direct held-input
comparison. Another full reconstruction or a higher iteration cap is not yet
supported by this evidence. The original container identity remains unknown.

Local evidence: `build/adam-feedback-2026-10-30/triage.json`,
`review-preservation-qa.json`, `tester-tools-tests.log`, and paired review copies.
Distribution: `build/pymacruise-review-tools-2026-10-30.zip`.
