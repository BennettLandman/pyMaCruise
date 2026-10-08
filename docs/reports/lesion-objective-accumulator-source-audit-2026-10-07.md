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
useful fixed-input comparisons for the all-positive-shape case. Its class
energy maps, auxiliary energy maps, running objective, and total objective do
not establish selected-source equivalence. The focused test and fixture README
now say so explicitly; plan row P39zzzzzzah is reopened.

## Python change

`run_lesion_toads_objective_membership_pass(..., source_faithful=True)` now
preserves the source's per-class and auxiliary energy values for active
classes and carries them through a zero-shape skip, which sets that class
membership to zero. The existing safer positive-shape
behavior remains available as the default. The current synthetic regression
checks current-versus-stale energy values and membership clearing, but it
is a source-contract regression, **not** a direct Java composite comparison.

Focused command:

```sh
./.venv/bin/python -m pytest -q \
  tests/test_lesion_objective_membership_pass.py \
  tests/test_lesion_descriptor_membership.py
```

Result on implementation commit `1843e87`: **6 passed in 4.63s**. No full
suite was run for this localized option. The Java method's complete
accumulator behavior remains unverified by a standalone composite harness;
that harness, source-faithful objective comparison, and selected production
integration remain open.

## Next validation

Build a small Java reference that preserves the original method-local energy
values and zero-shape branch, then compare every membership map, effective
energy used by the objective, denominator, running objective, and final total
against the Python source-faithful option. Include at least two dependent
interior voxels so both active reassignment and stale-value behavior are
observable. Do not reuse the previous per-voxel-energy JSON as the objective
golden.
