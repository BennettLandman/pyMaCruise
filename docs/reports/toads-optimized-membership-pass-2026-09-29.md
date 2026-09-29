# P18ar: composed optimized-tissue membership pass

## Scope

This adapter composes the previously validated optimized CSF/GM/WM stages
into one Python-native membership pass: image-wide inverse-variance weights,
assigned-class standard deviations, centroid-similarity coefficients,
optional PSF ridge bias, optimized `PowerTable` likelihood terms, atlas
inverse-shape-prior energy, six-face smoothing, relation factors, reciprocal
conversion, normalization, and source-ordered objective accumulation.

All membership classes in this helper are optimized tissue classes. Inputs
use the three optimized tissue channels from existing synthetic fixtures;
closest-label references to the additional fourth shape fixture class are
mapped to no match for this three-class slice. The helper runs one pass only:
centroid re-estimation, topology admissibility, atlas deformation, other
object types, iterative convergence, and complete segmentation remain
outside scope.

## Java source-formula references

The harness reuses committed 7×6×5×3 RAS 1 mm image and field volumes,
membership channels, categorical assignments, and transformed shapes from
earlier synthetic fixtures. It calls the actual selected Java `PowerTable`
class and combines the relevant `MultipleToadSegmentation` formulas. Four
committed JSON references cover field on/off crossed with PSF ridge on/off.
The selected source SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.

Regenerate the reference outputs with
`scripts/generate_java_toads_optimized_membership_pass_reference.sh`.

## Validation

`tests/test_toads_optimized_membership_pass.py` compares all stage
intermediates and final outputs. Image variance weights, global and
class-specific deviations, centroid similarities, and PSF biases match
exactly as float32. Data energy differs from Java by at most `2.98e-8` over
the four cases. Final memberships differ by at most `8.94e-8`; cumulative
objective differs by at most `1.20e-7`. The tests also verify untouched
source borders and reject invalid closest-label slot counts. The focused
module passes **5 cases**.

The full configured suite passes **348 tests in 62.10 seconds** on
implementation commit
`a7365b92891f6a7dce1a01e03d709e2e89be6097`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This is a composed
single-pass source comparison, not a full TOADS segmentation result.
