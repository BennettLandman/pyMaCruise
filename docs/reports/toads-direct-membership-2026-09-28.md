# TOADS direct reciprocal-energy membership update

**Date:** 2026-09-28  
**Roadmap:** P18ah complete  
**Validated implementation/fixture commit:** `e652c290efcf0c6c3d685f7282dd27e852c7a524`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_direct_membership_update` implements the
selected direct-relation update for optimized tissue classes after their
energy terms are assembled. At each interior voxel, a class energy greater
than `ZERO=1e-10f` is replaced by its float32 reciprocal; a value at or below
the threshold is replaced by `INF=1e10f`. The class reciprocals are added in
class order using float32, then divided by their sum when the denominator is
positive. The source objective uses `objectiveEnergy * membership *
membership`; the helper returns those per-class contributions. Outer-volume
memberships remain as initially supplied.

The selected `MultipleToadSegmentation.java` source SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
`ToadsDirectMembershipReferenceHarness.java` contains a formula excerpt for
the direct reciprocal, denominator, normalization, and objective operations.

## Synthetic fixture and results

The committed `tests/data/synthetic/toads_direct_membership/` directory
contains data-energy, objective-energy, and initial-membership float32
volumes, each 7×6×5×4 in 1 mm isotropic RAS `.nii.gz`. The energy phantom
includes zero, exact threshold, one float32 step below, and one float32 step
above `ZERO`. All values are verified against the generator.

Updated memberships and all per-class objective contributions match the Java
formula reference exactly after float32 conversion. The four threshold
locations choose the expected class, border memberships remain unchanged, and
border objective contributions stay zero. Six focused tests pass, including
invalid grid/threshold checks.

## Scope limits

This validates direct optimized-class inversion and normalization only.
P18ag separately validates the smoothing increment, but this slice does not
compose the asynchronous source membership sweep. Atlas shape priors, relation
factors, lesion/outlier mixtures, topology gates, and full TOADS evolution
remain open. The reference is a source-formula excerpt, not a full Java
segmentation output.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **297 passed in 51.50 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `e652c290efcf0c6c3d685f7282dd27e852c7a524`
- Run date: 2026-09-28
