# Full pytest report — P18 and plus compositions — 2026-10-01

## Result

- Implementation commit: `8d638e64efa45075e78c9e173d73e1ea86e0e129`
- Command: `./.venv/bin/pytest -q`
- Result: **993 passed, 0 failed, 0 skipped**
- Expected RAS-orientation warnings: **11**
- Elapsed time: **725.40 seconds** (12 minutes, 5 seconds)
- Data scope: committed synthetic fixtures and bundled fixed atlas resources only; no subject imaging data.

## Warning categories

The 11 warnings are the suite's expected notices for inputs normalized to RAS: two image-calculator inputs, one N4 input, two Dura-topology inputs, and six workflow-preflight inputs. The new plus-composition orientation test captures and asserts its LAS-to-RAS warning, so it does not add an uncaptured warning to this summary.

## Added coverage in this run

This run includes the P18cyu one-channel membership → selected registration/control → growth/relation/terminal branch composition and P18cyv continuation into the selected next thinning pass. The next-thinning Java reference consumes the Python-produced grow/relation state and independently compares the downstream queue and thinning outputs; it does not independently recompute growth or relations, and the selected order harness is metadata rather than execution of the Java wrapper.

P39 coverage now includes the selected post-loop output dispatch and the bridge from its `dura removal inputs` product through three transferred channels and shared crop into eight plus-ready RAS images. The Java comparison is limited to the identity RAS-grid path. The selected JIST orientation-only transfer does not generally equal the physical-affine transfer helper for differing origins, dimensions, or spacings.

P40 coverage includes the selected boundary → tissue fusion → ACE → GVF composition, with GVF receiving ACE's enhanced GM output, followed by the nested TGDM adapter on an explicit topology-corrected initial membership. Inputs and outputs are checked on one crop-ready RAS grid, including a non-RAS input warning/reorientation case. These are synthetic composition checks; topology correction and crop construction are supplied inputs, and complete MaCRUISEplus equivalence is not established.

The suite also covers P18cyt registration/control with final Java lambda and all class-bound comparisons. P59 adds scoped registration literature and code annotations; it identifies ANTsPy as an Aladin substitute and makes no numerical-equivalence claim.

## Scope and limitations

These additions remain bounded synthetic components and compositions. They do **not** establish full TOADS, full LesionTOADS, full MaCRUISEplus, or P37 whole-workflow numerical equivalence. P18 and P39 remain Started, P34 remains Blocked on complete TOADS/Dura behavior, and P41 plus-workflow assembly remains open. See the [P35 synthetic readiness inventory](p35-synthetic-readiness-inventory-2026-10-01-latest.md) and [plan](../../plan.md).
