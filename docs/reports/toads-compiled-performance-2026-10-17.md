# TOADS compiled execution and performance — 2026-10-17

## Change and numerical contract

The production exact-membership TOADS path now compiles three remaining expensive operations:

- Growing initialization and the complete min-heap loop, including 18/6 topology-LUT checks. The event-recording Python path remains available. Callback reporting runs between bounded batches without changing event order.
- Relation diffusion, closest-label diffusion, and boundary-distance seeding. Double buffers preserve synchronous reads, fixed face order, first/last candidate rules, iteration counts, borders, and carried byte values. The selected boundary-distance diffusion is a mathematical no-op and is not replaced with a different distance transform.
- Degree 1–3 exact Chebyshev bias-field accumulation, evaluation, and centering. Float32 intermediate rounding, source scan/reduction order, the float32-derived PowerTable step, and the existing JAMA LU solver remain unchanged. Strided multi-channel fields are still updated in place.

Numba is already a required dependency. No parallel reductions, relaxed floating-point operations, changed iteration schedules, or changed stopping settings are introduced. Each accelerated public component accepts `backend="python"` for the retained reference implementation. Growing uses the compiled path only when event snapshots are disabled. The zero-threshold NaN-priority behavior is preserved.

Implementation commits: `f2158a79` (kernels), `6c419308` (exact regression tests), and `6c3ad57a` (NaN edge behavior and benchmark driver). Later report/driver changes do not change the numerical kernels.

## Controlled synthetic performance

The committed 36×45×34 folded-atlas phantom ran through its natural stop at iteration 19 with the compact atlas, padding 2, one registration level, and one loop registration iteration. Each backend ran in a fresh process on the same host with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Compilation/cache loading and output I/O are included in total times. Other validation jobs were running; these are observed wall times rather than isolated CPU benchmarks.

| Work | Python reference | Compiled | Speedup |
|---|---:|---:|---:|
| Complete run, including alignment/output | 36.464 s | 8.287 s | 4.40× |
| Growing, 19 calls | 13.377 s | 0.205 s | 65.1× |
| Relation map, 19 calls | 2.311 s | 0.112 s | 20.6× |
| Closest labels, 19 calls | 2.032 s | 0.079 s | 25.8× |
| Boundary distances, 19 calls | 0.230 s | 0.012 s | 19.2× |
| Degree-three field fits, 8 calls | 9.310 s | 0.022 s | 423× |

The classification files are byte-identical (SHA-256 `652ac6c2bbba0faeedecf56d033ee9ec6a9f668f889f13d0e273c3ace55aedd1`) and match the fresh pre-change baseline exactly. A separate uninstrumented loop comparison measured 31.169 s before and 5.305 s after (5.88×), excluding alignment and output. Those numbers describe a different timing boundary from the full-run row.

This compact benchmark uses different atlas/padding settings from the committed nineteen-iteration Java fixture. It is a matched Python execution comparison. The existing full-atlas Java membership/pose/label regression supplies the separate source check; the compact benchmark must not be compared directly with that fixture's classification.

## Human replay

A fresh first-TOADS replay used the saved case-100260 command from `build/sample-singularity-base-source-profile-2026-10-15/run.json`: the exact saved unmasked-N4 input, brain mask, atlas descriptor, four-level registration settings, 50-iteration cap, and unchanged stopping thresholds. Outputs went to a new ignored directory.

| Timing boundary | Saved Python baseline | Fresh optimized run | Observed speedup |
|---|---:|---:|---:|
| Load inputs through loop completion | 3,012.1 s (50.20 min) | 809.9 s (13.50 min) | 3.72× |
| Segmentation loop, excluding alignment | 2,675.8 s (44.60 min) | 497.5 s (8.29 min) | 5.38× |

The complete fresh CLI invocation, including output writing, took **810.449 s**. The baseline elapsed values come from the saved same-settings run, not a new contemporaneous unoptimized human replay. The controlled synthetic A/B above provides the fresh reference-versus-compiled comparison.

Both human runs stop naturally at **iteration 20**, with terminal energy ratio **4.29838547e-5**. The printed energy ratios match at all 20 iterations. Classification differs at **zero voxels**, and all five emitted Dura-input files (classification, CSF, GM, WM-fill, WM-mask) are **byte-identical**, with identical data types and physical grids. The saved baseline had already matched selected Java labels on these fixed inputs; this fresh optimization run does not repeat the full external Java/MATLAB workflow.

The first growing pass exhausts the same **1,548,285 events**. Its queue-processing portion is approximately 0.7 s versus the saved 29.4 s, while initialization/cache loading is additional. Its relation refresh is approximately 1.0 s versus the saved 30.6 s. Across the full fresh human command, growing totals 13.090 s over 20 calls; relation-map, boundary-distance, and closest-label operations total 20.139 s; all 14 field fits total 1.964 s, including 9 degree-three fits in 1.446 s. Repeated rigid registration, initialization, and membership/statistics now dominate the remaining runtime.

Process high-water RSS was **7.207 GB**; the one-second watchdog sampled a peak of **6.702 GB**, below its 12 GB limit. The higher process high-water value illustrates that polling can miss brief peaks. Machine-readable timings and output comparisons are in [the evaluation JSON](toads-compiled-performance-2026-10-17.json).

## Verification

Focused checks passed:

- Growing terminal state versus selected Java heap/LUT fixtures and degree-two/three direct Java field references: 12 passed in 7.69 s.
- Selected Java relation/closest-label/boundary-distance fixtures: 30 passed in 2.53 s.
- Additional exact Python-versus-compiled checks, including degree-one/two/three intermediate matrices and coefficients, strided fields, synchronous diffusion, borders, signed-byte carried state, tied priorities, input nonmutation, progress batching, and zero-threshold NaNs: 18 passed in 4.38 s.

The full regression suite passed **1,670 tests, 4 skipped, 0 failed**, with 21 expected orientation warnings in **721.40 s**. Peak sampled RSS was **8.475 GB** under a 16 GB watchdog. The extra NaN edge test was added after full-suite collection and passed in the focused 18-test run above. See [full-suite evidence](full-pytest-2026-10-17-toads-compiled.md) for the command, collection boundary, and corrected import-path invocation.

## Reproduction

From the repository root, run each command separately:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/benchmark_toads_optimizations.py --reference \
  --output-dir build/toads-performance-reference

PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/benchmark_toads_optimizations.py \
  --output-dir build/toads-performance-compiled \
  --compare build/toads-performance-reference/classification.nii.gz
```

For the human stage, replace the compact invocation with:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -u scripts/benchmark_toads_optimizations.py \
  --source-run build/sample-singularity-base-source-profile-2026-10-15/run.json \
  --output-dir build/toads-performance-human
```

The driver refuses to overwrite prior results, redirects all stage outputs, records component times and classification hashes, and compares human labels with the saved same-settings classification. `--stage toads-second` selects the other recorded pass. Raw logs and measurements for this evaluation remain under ignored `build/toads-performance-2026-10-17/`; no human images are committed.

This is an execution optimization and selected-stage regression, not new whole-workflow source equivalence or anatomical acceptance. Source-method coverage estimates remain unchanged.
