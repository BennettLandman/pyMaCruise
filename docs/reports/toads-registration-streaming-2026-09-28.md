# TOADS rigid-registration streaming slice

**Date:** 2026-09-28  
**Status:** P18v complete; source-level registration equivalence remains open.

## Streaming versus dense reference

`stream_toads_registration_energy` and
`stream_toads_registration_coefficients` in
`src/pymacruise/toads_registration_sampling.py` sample each required atlas
class and immediately accumulate its contribution. The optimizer uses these
functions instead of materializing the complete prior, derivative, and
Jacobian fields. The dense `sample_toads_registration_fields` helper remains
available for small-grid comparisons.

On a committed anisotropic 4×3×3 synthetic input grid with three classes, a
nonzero rigid transform, different image/template resolutions, and distinct
class bounds, streaming and dense paths produced exactly equal float32 energy,
cost, normalization, gradient, diagonal Hessian, norm-gradient, and norm-Hessian
for both full-grid and offset/subsample scans. The same test excludes one
registered class. The optimizer determinism and stop-condition tests now
execute through the streamed path. Tests: `tests/test_toads_registration_sampling.py`.

## Synthetic time and scratch-memory probe

Reproduce from the repository root with:

```sh
./.venv/bin/python scripts/benchmark_toads_registration_streaming.py --sizes 6 8 10
```

Each measurement includes one coefficient scan plus one energy scan over a
three-class deterministic synthetic atlas. Input arrays are allocated before
timing and memory tracing. Times below are untraced; Python-visible temporary
allocation is measured in a separate call with `tracemalloc`.

| Grid | Voxels | Coefficient + energy time | Traced temporary peak |
|---:|---:|---:|---:|
| 6³ | 216 | 0.1053 s | 0.0058 MiB |
| 8³ | 512 | 0.2481 s | 0.0057 MiB |
| 10³ | 1,000 | 0.4903 s | 0.0059 MiB |

The near-constant scratch peak is consistent with voxel-at-a-time operation;
it excludes preallocated input atlas/membership arrays and is not a process
RSS measurement. Elapsed time grows roughly with sampled voxel count. At this
rate, direct Python sampling is still too slow for clinical-size volumes or
the repeated trial evaluations in rigid optimization. No real imaging data
was used. Vectorized source-order sampling/accumulation or another measured
implementation improvement is still required before production-size use.

## Evidence limits

Streaming is compared exactly with the Python dense composition; its transform,
interpolation, derivative, Jacobian, and coefficient components have separate
synthetic Java API/source-formula comparisons recorded in P18i/P18o–P18s. No
single full `MultipleToadDeformableAtlas` execution is available in the local
reference harness, so this is not full-Java registration equivalence, known-
transform alignment validation, or TOADS segmentation validation. Atlas
descriptor versus NIfTI spacing is also unresolved and no bundled atlas was
executed through this optimizer.

Full repository validation for the P18v implementation commit:

- Command: `./.venv/bin/pytest -q`
- Result: **252 passed in 70.82 seconds**
- Tested commit: `0fe1161d6d5cf8f4af61b13346801ce9236d49dc`
- Run date: 2026-09-28
