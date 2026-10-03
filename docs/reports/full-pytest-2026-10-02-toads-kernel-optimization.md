# Full synthetic suite: TOADS kernel optimization and cycle trace

- **Command:** `./.venv/bin/pytest -q`
- **Tested implementation commit:** `f3b874b`
- **Result:** 1,130 passed, 0 failed, 0 skipped; 11 expected RAS-orientation warnings
- **Elapsed time:** 1,213.79 seconds (20:13.79)

The suite includes bitwise scalar-oracle coverage for vectorized shape-relation
factors, selected-Java default-atlas membership-energy checks for batched
optimized data terms, the registered-loop selected-constructor setting
`do_growing=true`, a deterministic Python-only registered-state trace, and the
restartable synthetic base workflow. The existing two-update Java comparison
still checks the selected membership/statistics update seam. None of this is a
Java whole-outer-loop comparison or real-image validation.

## Synthetic TOADS diagnostics

Run the registered state trace with:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_state.py \
  --maximum-iterations 4 --output build/toads-python-state-trace-4.json
```

It writes JSON summaries and hashes at each iteration boundary. The report
identifies itself as Python-only synthetic evidence, with both source
equivalence and real-data validation set to false. For this 13³ fixture, the
membership state changes each cycle while the field and rigid transform stay
fixed. This narrows the investigation but does not identify a source mismatch.

The 50-cycle run used:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --output build/toads-loop-50-after-kernel-optimizations.nii.gz \
  --maximum-iterations 50 --convergence-threshold 1e-4
```

It took 62.496s and reached the iteration cap with terminal ratio `0.0385285541`;
the `1e-4` threshold was not met. The objective settled into alternating bands
near 64.85 and 62.40. The uint8 output remained on the native RAS grid and had
labels 1–11 (SHA-256
`eac4db1fc1a32b29f70e7d47eeea6bcfa1df020d34a4ffa9dbcaf799e65b5453`).

## Numerical-preserving performance work

- Ten optimized tissue-class data terms now use one ordered call; MASK remains
  separate. Against the synthetic selected-Java fixture, class-energy and
  complete membership outputs matched exactly. A generated 13³ direct term
  benchmark improved 10.2%; the whole tiny atlas pass improved only about 0.5%,
  so the ordered in-place membership sweep remains the main per-cycle cost.
- Shape-relation-factor loops now compute direct/neighbor masks over
  voxel/class arrays while retaining float32 operation order. A deterministic
  scalar oracle matched all factors and suppression flags bit-for-bit across
  clamp, zero, overlap, and signed-zero cases. On a generated 13³×11 array,
  the operation took about 0.0003s versus 0.0834s for the scalar oracle; this
  microbenchmark is not an end-to-end workflow speed claim.

These compose already represented Java operations; **the suite adds zero newly
counted selected Java source lines**. The planning estimate remains about
3,464–4,864 of 13,651 selected lines (25–36%); it is method-level scope
accounting, not executed source-line coverage or workflow completion.

## Follow-up after this suite

After the tested commit, the selected `prior_scale` default was corrected from
0.05 to the wrapper's atlas range 0.2. Focused registered-loop plus restartable
base tests passed **4 tests in 80.77s** on the implementation in the subsequent
project commit. The source loop was also compiled and run on the synthetic
fixture: it met its `1e-4` ratio at iteration 45, but its MIPAV working grid is
35×34×33, compared with Python's 19×17×15. Python with the corrected atlas
range still reached its 50-iteration cap at ratio 0.13936983. These are not
like-for-like convergence results; see the [VandyCRUISE runtime report](toads-vandycruise-synthetic-runtime-2026-10-02.md).

The next integration milestone is a repeatable source-side state trace through
preprocessing, startup, exact update, registration, growth, relation refresh,
and next thinning, with matching work grids and parameters. Full base-workflow
convergence, plus integration, and real-data validation remain open.
