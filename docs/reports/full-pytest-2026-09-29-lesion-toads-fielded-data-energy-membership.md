# Full test suite: LesionTOADS field-corrected data-energy membership integration

- Date: 2026-09-29
- Tested implementation commit: `499de268c1ca4c1625d46d8b840e29ad673017c7`
- Command: `.venv/bin/pytest -q`
- Result: **457 passed in 56.82 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes two selected-Java one-voxel compositions. Both multiply each
synthetic channel by its own nonunit intensity-field sample before computing
class and lesion residuals; the BLACKHOLE case also compares black-hole
residuals. The resulting terms pass through smoothing, shape/similarity
priors, direct/prior conversion, distance modulation, normalization, and
local-energy accumulation. All compared values match Java float32 bits.

The fielded fixture is committed under
`tests/data/synthetic/lesion_toads_fielded_data_energy_membership/`. A
separate no-argument harness invocation was also byte-compared against the
existing P39zk fixture to confirm the harness extension does not alter the
unfielded reference. Spatial stencils and shape/similarity vectors remain
fixed synthetic inputs; this is scalar integration, not a volume run.
