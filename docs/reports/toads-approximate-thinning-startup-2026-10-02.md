# Approximate-membership TOADS thinning startup

## Scope

Implemented the `useApprox` branch of VandyCRUISE `MultipleToadSegmentation.initThinningLabels` as `initialize_toads_approximate_thinning` in `src/pymacruise/toads_approximate_thinning_startup.py`. The method now supports the sparse top-k memberships produced by the approximate-membership path when preparing a thinning queue.

The port preserves the selected Java behavior: interior-only scan order; segmentation-to-skeleton copy; float32 accumulation of a current class's represented top-k memberships; initial availability reset inside the volume; six-face boundary detection; the source's 26-neighbor endpoint count; and negative reciprocal top-k membership priorities using its `ZERO` floor. The initialized max-heap uses the shared source-compatible heap implementation.

## Evidence

- Reference: a compact Java harness reproduces this selected method arm on a deterministic 5 x 5 x 5, three-class, two-slot synthetic grid. It records class-volume float bits, boundary counts, skeleton and availability arrays, and each queued priority and coordinate.
- Test: `tests/test_toads_approximate_thinning_startup.py` compares every recorded output to the Java harness reference and checks invalid repeated class slots.
- Focused validation: `./.venv/bin/pytest -q tests/test_toads_approximate_thinning_startup.py tests/test_toads_thinning_transition.py tests/test_toads_initial_thinning_startup.py` — **5 passed in 6.58 seconds**.
- Source file SHA-256: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Java harness SHA-256: `6369e2a3289e0295804f512c7728327345f37660a0078517b26b60a7897f796a`.

## Limits

This establishes a match for synthetic method-level startup state, not a full approximate TOADS iteration or segmentation equivalence. Inputs must already contain normalized sparse memberships and unique top-k class indices; this function does not generate them. Queue items are compared as a complete priority/payload multiset here; heap tie-removal order remains covered by the shared heap tests. The fixture is not an imaging-data validation.
