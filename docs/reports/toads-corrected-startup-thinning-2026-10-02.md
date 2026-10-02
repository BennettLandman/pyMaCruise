# Field-corrected registered startup through first TOADS thinning (P18cyax)

## Ported behavior

`compute_registered_toads_corrected_initial_membership` now continues the
registered startup beyond its degree-three correction fit and field-aware
statistics/membership pass. It passes the updated memberships, the registered
class adjacency, and transformed shape priors into the selected first
`initThinningLabels`/`propagateThinning` transition. The relation state is
reused from registration startup: the wrapper refreshes relations after
attaching the field, but no memberships or labels change before that refresh.
Volume ratio and thinning limit remain explicit settings with selected
defaults of `0.5` and `3.0`.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_field_refit.py tests/test_toads_initial_statistics_thinning.py tests/test_toads_default_continue_thinning.py tests/test_toads_thinning_transition.py tests/test_toads_initial_relation_state.py tests/test_toads_default_first_thinning.py
```

Result: **13 passed, 0 failed, 0 skipped in 31.57 seconds** on implementation
and test commit `7d83dc95e90762baa1a334ae8496e85c0202f25e`.

The integration test compares segmentation, skeleton, and processed-event
count with direct calls using the same field-corrected membership, relation,
and shape arrays. The field fit has its selected Java formula comparison;
relation, membership, and thinning primitives have their own source-reference
tests. This composition is not an independent combined Java wrapper oracle,
does not include the following growing pass, and is not full TOADS equivalence.
The full repository suite has not yet been rerun, so P18cyax remains Started.
