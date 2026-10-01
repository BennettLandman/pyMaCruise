# Full pytest report — source-matched continuations — 2026-10-01

## Result

- Implementation commit: `c747e5b0e7ae4f9423d0ee771ee39bd66ebf4780`
- Command: `./.venv/bin/pytest -q`
- Result: **986 passed, 0 failed, 0 skipped**
- Expected RAS-orientation warnings: **11**
- Elapsed time: **690.19 seconds** (11 minutes, 30 seconds)
- Data scope: committed synthetic fixtures and bundled fixed atlas resources only; no subject imaging data.

## Warning categories

The 11 warnings are the suite's expected notices for inputs that are reoriented to RAS: two image-calculator inputs, one N4 input, two Dura-topology inputs, and six workflow-preflight inputs. No warning indicated a failed test.

## Scope

This full run includes the source-matched TOADS and LesionTOADS component continuations at this commit, including the Dura attached/no-dura wrapper-stage reference, the one-channel TOADS degree-one field/statistics/membership chain, selected next-loop thinning, LesionTOADS output dispatch, the Dura-input terminal branch, and the NaN-gated thinning regression. These remain bounded synthetic component or wrapper-sequence references. They do **not** demonstrate complete TOADS or LesionTOADS workflow equivalence, satisfy the P37 whole-workflow numerical gate, or unblock P34.

The P35 synthetic readiness inventory was refreshed against this same implementation commit in [the inventory report](p35-synthetic-readiness-inventory-2026-10-01.md). The subsequent documentation-only commit records both reports and the corresponding roadmap evidence.
