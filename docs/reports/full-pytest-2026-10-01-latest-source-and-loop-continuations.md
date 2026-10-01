# Full pytest report — latest source and loop continuations — 2026-10-01

## Result

- Implementation commit: `5cbd3856129fa22fca524e3bc75bf3eb8f56eadd`
- Command: `./.venv/bin/pytest -q`
- Result: **1,022 passed, 0 failed, 0 skipped**
- Expected RAS-orientation warnings: **11**
- Elapsed time: **1,330.98 seconds** (22 minutes, 10.98 seconds)
- Data scope: committed synthetic fixtures and bundled fixed atlas resources only; the human MRI run was not used.

## Warning categories

The 11 warnings are the expected RAS reorientation notices: two image-calculator inputs, one N4 input, two Dura-topology inputs, and six workflow-preflight inputs. No warning indicated a failed test.

## Scope and limits

This run includes the latest bounded P18 and P39 Java-reference continuations through selected TOADS registration/field/membership/loop/output paths, LesionTOADS later membership and registration/control paths, terminal label routing, and terminal crop/buffer materialization. It also exercises the source-order-specific anisotropic TOADS registration coordinate path. These are selected synthetic operations and wrapper slices; they do **not** establish complete TOADS or LesionTOADS workflow equivalence, satisfy whole-workflow tolerance gates, or unblock P34. No subject imaging data were used.

The complete suite is reproducible at the implementation commit with the command above. The following documentation-only commit updates the plan and readiness inventory to point to this report.
