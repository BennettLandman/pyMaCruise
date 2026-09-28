# TOADS registration vectorization

**Date:** 2026-09-28  
**Roadmap:** P18w complete  
**Tested implementation commit:** `dad38b855a480638eb717c0d90421b5c8aa9de38`

## Numerical validation

`vectorized_toads_registration_energy` and
`vectorized_toads_registration_coefficients` use NumPy operations over one
input x-slab at a time. Transform mapping, closest-edge trilinear sampling,
trilinear derivatives, and the rigid Jacobian are vectorized over each slab.
Coefficient terms are laid out in source x/y/z/class order, then accumulated
with float32 cumulative sums seeded by the preceding slab's totals.

Two synthetic anisotropic grids cover distinct image/template spacings,
nonzero rigid rotations and translations, class-specific support bounds, and
template-edge/out-of-bounds behavior. For full scans and offset/subsample
scans, the vectorized energy, cost, norm, six gradients, six diagonal Hessians,
six norm-gradients, and six norm-Hessians match both the scalar streaming path
and dense field reference exactly as float32 values. The integrated optimizer
uses the vectorized kernels; deterministic history, accepted updates, and
outer stop conditions pass on synthetic fixtures. No clinical data are used.

## Synthetic benchmark

Reproduce from the project root:

```sh
./.venv/bin/python scripts/benchmark_toads_registration_streaming.py --sizes 6 10 16
```

Each time combines one rigid coefficient scan and one objective scan for three
synthetic classes. Input arrays are allocated before timing. Scratch memory is
measured separately with `tracemalloc` during one vectorized coefficient scan;
this is Python-visible temporary allocation, not process RSS, and excludes the
input arrays.

| Grid | Voxels | Scalar streaming | Vectorized | Speedup | Traced scratch peak |
|---:|---:|---:|---:|---:|---:|
| 6³ | 216 | 0.10364 s | 0.00689 s | 15.05× | 0.0323 MiB |
| 10³ | 1,000 | 0.48233 s | 0.01176 s | 41.02× | 0.0699 MiB |
| 16³ | 4,096 | 2.03102 s | 0.02076 s | 97.83× | 0.1580 MiB |

The fixed-size x-slab reduces Python-loop overhead while preserving measured
source-order sums exactly. These small-grid times do not prove runtime on the
bundled atlas or a full TOADS run. The optimizer still repeats sampling for
energy trials, atlas-grid resolution disagrees between the descriptor and
image metadata, and no full Java registration comparison has been run.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **253 passed in 75.48 seconds**
- Python: 3.13.2; pytest: 9.1.1
- Tested commit: `dad38b855a480638eb717c0d90421b5c8aa9de38`
- Run date: 2026-09-28
