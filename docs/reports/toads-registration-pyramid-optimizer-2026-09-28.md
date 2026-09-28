# Coarse-to-fine TOADS rigid registration

**Date:** 2026-09-28  
**Roadmap:** P18z complete  
**Implementation commit:** `23584b675c57b72b6014930bc28a687de63159ee`

## Source behavior ported

The supplied VandyCRUISE
`MultipleToadDeformableAtlas.registerShapesExactPyramid` starts at
`2**(levels - 1)`, runs each coarser scale down to 2, and omits scale 1. At
each level it block-mean subsamples memberships and each class prior directly
from its full-resolution array. It resets the line-search lambda to its
initial value, runs the exact-membership optimizer, and carries the updated
rigid transform into the next finer level. The source rigid coordinate and
Jacobian methods keep the original full-grid centers and divide mapped
template coordinates and parameter derivatives by the current scale.

`register_toads_atlas_rigid_pyramid` in
`src/pymacruise/toads_registration_sampling.py` wires together the previously
validated P18x/P18y kernels and the existing rigid optimizer. Each level is
built from the original membership and atlas arrays. It passes explicit
full-resolution centers, per-scale parameters, and full coarse-grid class
bounds into the optimizer. This follows the selected exact source registration
loop, which visits each coarse image voxel and each registered class, and does
not use the separate fixed-transform support bounding boxes. Scale 1 returns
an empty level history and leaves the initial transform unchanged.

## Synthetic validation

The integration test reads committed float32 4D RAS NIfTI fixtures in
`tests/data/synthetic/toads_registration_pyramid_optimizer/`: 17×15×13×3
memberships and 19×17×15×3 atlas priors. It verifies their voxel values
against NumPy `default_rng(2)` and confirms positive RAS axis codes. With
nonidentity initial rigid parameters and anisotropic resolutions, it runs scales 4 then 2,
accepts updates at the first level, passes that exact resulting transform as
the next level's initial state, and reproduces both complete level histories
and the final transform exactly on repeat. The P18x Java formula fixture
independently verifies the block means, and the P18y Java formula fixture
checks scaled mapping and Jacobian values at scales 1, 2, and 4.

A separate scale-2 fixture compares the vectorized optimizer kernels with
per-voxel scalar shape sampling, prior derivatives, rigid Jacobians, and the
existing exact-membership accumulation kernel. Energy is exactly equal. Raw
coefficient arrays pass `rtol=2e-7`, `atol=3e-10`; maximum absolute difference
is `1.49e-8`, attributable to float32 rounding in the vectorized accumulated
coefficients. This small coefficient deviation is recorded rather than hidden
inside an alignment tolerance.

The tests also verify the source's one-level/no-scale-1 behavior. The fixtures
are synthetic probability fields, not clinical scans; bundled atlas files are
not read by this test.

## Scope limits

P18z validates only the Python optimizer's multilevel wiring and its component
kernels. It does not execute the private Java registration method, establish
atlas alignment accuracy, run deformable TOADS tissue segmentation, resolve
the bundled atlas spacing discrepancy, or unblock end-to-end P34. Those gaps
remain tracked in [the TOADS registration notes](../toads-atlas-registration.md)
and [the main plan](../../plan.md).

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **266 passed in 56.94 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested implementation commit: `23584b675c57b72b6014930bc28a687de63159ee`
- Run date: 2026-09-28
