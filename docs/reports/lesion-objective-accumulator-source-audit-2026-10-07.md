# LesionTOADS objective accumulator source audit

## Finding

The prior P39zzzzzzah harness overstated its scope. In the selected VandyCRUISE
`LesionToadSegmentation.computeMemberships()` implementation,
`energy` is a method-local `float[]` initialized once at line 994, while
`lesion_energy` and `blackHole_energy` are initialized at line 995, then
updated as `energy[k] += ...` while the x/y/z scan runs (for example, line
1363). The objective later uses that accumulated value at each voxel (line
1547). The class loop also gates updates on `shape[k] > 0` (line 1084); when
the gate fails, the later `else` writes zero class membership, while the energy
accumulator is left unchanged. Auxiliary `lesion_energy` and
`blackHole_energy` are also method-local accumulators (lines 1000–1001).

The earlier harness instead stored a fresh per-voxel energy array and reduced
that array into its objective. Its membership and denominator values remain
useful fixed-input comparisons for the all-positive-shape case. Its class
energy maps, auxiliary energy maps, running objective, and total objective do
not establish selected-source equivalence. The focused test and fixture README
now say so explicitly; plan row P39zzzzzzah is reopened.

## Python change

`run_lesion_toads_objective_membership_pass(..., source_faithful=True)` now
preserves the source's per-class and auxiliary energy accumulation across the
full scan, and its zero-shape skip sets that class membership to zero while
leaving the energy accumulator unchanged. The existing safer positive-shape
behavior remains available as the default. The current synthetic regression
checks accumulated-versus-stale energy values and membership clearing, but it
is a source-contract regression, **not** a direct Java composite comparison.

Focused command:

```sh
./.venv/bin/python -m pytest -q \
  tests/test_lesion_objective_membership_pass.py \
  tests/test_lesion_descriptor_membership.py
```

Result on implementation commit `1f411b6`: **6 passed in 0.65s**. No full
suite was run for this localized option. The Java method's complete
accumulator behavior remains unverified by a standalone composite harness;
that harness, source-faithful objective comparison, and selected production
integration remain open.

## Next validation

Build a small Java reference that preserves the original method-local
accumulators and the zero-shape branch, then compare every membership map,
effective energy used by the objective, denominator, running objective, and
final total against the Python source-faithful option. Include at least two
dependent interior voxels so both accumulated and stale-value behavior are
observable. Do not reuse the previous per-voxel-energy JSON as the objective
golden.
