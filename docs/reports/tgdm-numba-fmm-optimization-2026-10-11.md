# TGDM distance reinitialization optimization — 2026-10-11

## Change

TGDM signed-distance reinitialization now uses a dedicated compiled Numba
ordered-fast-marching kernel in `src/pymacruise/tgdm_fmm.py`. It retains the
JIST source scan order, six-neighbor update order, heap tie behavior,
subvoxel seed arithmetic, float32 queue priorities, and float32 output storage.
The weighted ACE distance/characteristic path remains on its existing Python
implementation.

This targets the central FMM that prevented the previous human TGDM-only stage
replay from completing within its 15-minute cap. This report establishes
kernel-level equivalence and runtime only; it is not a human TGDM result.

## Validation

- `./.venv/bin/pytest -q tests/test_tgdm_fmm.py tests/test_tgdm.py`: **29 passed
  in 4.23 s** on the current working tree.
- The implementation author also ran
  `tests/test_tgdm.py tests/test_tgdm_fmm.py tests/test_toads_cruise_core.py`:
  **32 passed in 8.04 s**; focused ACE regressions plus FMM tests: **5 passed
  in 1.78 s**; original-JIST JAR tests: **4 passed in 0.59 s**.
- Four saved original-JIST `DistanceField` fields match exactly. The compiled
  path also matched the prior Python path bit-for-bit on a folded 48×47×43
  field and 20 seeded random boundary/sign fields; finite-band output equals
  clipping the full field.
- On a synthetic 161×195×154 folded field, the compiled call was **1.395 s
  cold** and **0.130 s warm**, compared with **19.07 s** for the prior Python
  path. The compiled call's measured process peak RSS was **380.7 MiB**. No
  OS-level memory limit was enforced in that benchmark.

These are source-reference and synthetic performance checks. Human integration
must still be rerun after the optimization.

## Next validation

Rerun the previously bounded human TGDM stage-isolation script with a 12 GiB
RSS watchdog. If it completes and returns stage fields, compare the saved
Python stages with the prepared source-selected Java TGDM harness on the exact
same cropped inputs and parameters. Preserve the output under ignored
`build/`; do not commit human images or derived images.
