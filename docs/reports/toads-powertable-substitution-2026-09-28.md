# TOADS PowerTable substitution comparison

**Date:** 2026-09-28  
**Roadmap:** P18ad complete  
**Validated implementation/fixture commit:** `b0dd1fe4aafd9f78b66b137fa0ed04c5630c65fb`

## Selected source contract

The selected `MultipleToadSegmentation` builds three `PowerTable` objects for
the optimized CSF, GM, and WM likelihood exponents. For the selected T1 MPRAGE
descriptor, those exponents are `[0.8, 1.0, 1.0]`. Its constructor call passes
`0.000001f` as a Java float which is then promoted to the `double` field:
`9.999999974752427e-7`. For a cached exponent and in-range input, the actual
`PowerTable.lookup` computes
`(int)(((p - pmin) / pstep) + 0.5)` and evaluates the cached power at that
rounded table sample. Out-of-range values and cache misses call `Math.pow`
directly.

The complete selected source class is preserved with its author/version
comments at
`tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java`.
Its SHA-256 is
`2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`.
`ToadsPowerTableReferenceHarness` compiles that actual class and writes the
committed `java_reference.json` fixture.

## Results

The 1,042 deterministic cases cover endpoints, tiny values, values around
lookup half-step boundaries, and 1,024 `Random(1804)` values across `[0,4]`.
Each case includes a three-class likelihood vector and the source float32
inverse-energy normalization.

- Direct unquantized `numpy.power` maximum absolute power error: `8.96e-6`.
- Maximum relative power error: `0.426`, concentrated at very small inputs.
- For input powers `p >= 0.01`, maximum absolute power error: `8.07e-7`.
- Direct NumPy maximum normalized-membership error: `0.4847`.
- Direct NumPy class-winner disagreements: `3 / 1,042`.
- `toads_powertable_lookup` rounds inputs using the Java float32-derived step,
  then calls NumPy `power`; normalized memberships differ by less than `1e-6`
  on all cases, and no winning labels change.

The scalar power values can look close while the normalized inverse-energy
membership changes sharply near the source table's zero boundary. Therefore
direct unquantized NumPy power is not used for this source path. The small
wrapper in `src/pymacruise/toads_math.py` reuses NumPy's power operation while
retaining Java's index rounding and direct-power fallback behavior.

## Scope limits

This comparison validates the table lookup and one source-formula
three-class inverse-energy normalization. It does not run `computeExactMemberships`,
the complete optimized tissue likelihood with smoothing/atlas priors, fast
marching, topology-constrained updates, or full TOADS segmentation. Those
remain open P18 work. No clinical data are used.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **270 passed in 52.01 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested commit: `b0dd1fe4aafd9f78b66b137fa0ed04c5630c65fb`
- Run date: 2026-09-28
