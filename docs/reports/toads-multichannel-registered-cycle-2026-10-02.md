# Multi-channel TOADS continuation through the next thinning pass (P18cybc)

## Ported behavior

The multi-channel exact-membership continuation now feeds the selected
wrapper sequence after memberships: compute the energy ratio, run rigid atlas
registration, apply convergence and iteration-limit checks, grow when the
branch continues, refresh relations/boundary distances/closest labels, and
start the next thinning pass using the supplied transformed shape priors.
Existing selected registration, control, growth, relation, and thinning
kernels are reused.

## Validation

Focused test:

```text
./.venv/bin/pytest -q tests/test_toads_multichannel_carried_degree_iteration.py::test_multichannel_memberships_continue_through_registration_growth_and_next_thinning
```

Result: **1 passed, 0 failed, 0 skipped in 11.90 seconds**. A subsequent run
of the full affected module passed **3 tests in 25.55 seconds**.

The case supplies committed synthetic 21³ RAS state, exercises a real
one-channel atlas registration over two-channel memberships, and continues
through a non-stopping growth/relation branch into the next thinning. It
checks membership object identity at registration, presence of registration,
continuing control, growth and refreshed relations, stage order, and equality
of the next thinning result with a direct selected component call. This is a
composition check over separately source-checked operations; it is not a
combined Java wrapper execution or full TOADS equivalence.

## Scope limits

The caller still supplies transformed shapes, geometry, prior loop energy,
hard-label/skeleton/ordering/volume state, and atlas arrays. Stop-time output
transfer, initialization, field/statistics estimation for the next complete
iteration, and full workflow equivalence remain separate.
