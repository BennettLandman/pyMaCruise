# LesionTOADS objective accumulator source audit

## Finding

The prior P39zzzzzzah harness overstated its scope. In the selected VandyCRUISE
`LesionToadSegmentation.computeMemberships()` implementation,
`energy` is a method-local `float[]` initialized once at line 994; the two
auxiliary energy values are initialized at line 995. For each active class,
the data-term branch assigns `energy[k] = num` (lines 1252–1261), then
same-voxel smoothing/prior terms add to it (including line 1363). The class
loop is gated by `shape[k] > 0` (line 1084), so a zero-shape class skips both
the assignment and membership calculation; the later `else` writes zero class
membership. That skipped class therefore retains the previous energy value
for the objective at line 1547. Active lesion/black-hole classes similarly
assign their auxiliary data energies, then add same-voxel spatial terms.

The earlier harness instead stored per-voxel energy arrays and reduced those
arrays into its objective. Its membership and denominator values remain
useful fixed-input comparisons, but its class/auxiliary energies and objective
did not establish selected-source equivalence. The fixture now includes
zero-shape skips for all four supported class types. The corrected Java
excerpt harness retains the prior method-local energy value on a skip, resets
active-class energy from the voxel's data term, adds same-voxel spatial terms,
and compares every returned map and objective bit against Python.

## Python change

`run_lesion_toads_objective_membership_pass(..., source_faithful=True)` is now
the default and preserves the source's per-class and auxiliary energy values
for active classes and carries them through a zero-shape skip, which sets that
class membership to zero. The opt-out mode requires strictly positive interior
shapes. The current synthetic regression
is a direct selected-Java-excerpt composite comparison for the declared
four-class fixed-input slice.

Focused command:

```sh
./.venv/bin/python -m pytest -q \
  tests/test_lesion_objective_membership_pass.py \
  tests/test_lesion_descriptor_membership.py \
  tests/test_lesion_descriptor_iteration_slice.py \
  tests/test_lesion_continuing_iteration.py
```

Focused result on implementation commit `e774747`: **10 passed in 6.10s**. The pinned Java harness SHA-256 is `bbf086b2d69dc22d95454f43e4fd73d549538f6d8b081ceaa4910ae9794240ab`.
The full suite on the same commit passed **1,234 tests, 0 failed, 0 skipped,
11 expected RAS warnings in 1,675.72s**, with **10.831 GB** peak sampled RSS
under a 16 GB watchdog; see the [full-suite report](full-pytest-2026-10-07-lesion-objective-source-match.md).
Plan row P39zzzzzzah is closed for this bounded fixed-input slice. Full
LesionTOADS workflow integration and equivalence remain open.

## Next validation

Continue P39 by composing the validated pass with active statistics, atlas
sampling, registration, and later iterations; keep those integration claims
separate from this fixed-input result. Do not reuse the prior harness's
per-voxel-energy JSON as source-equivalence evidence for other profiles.
