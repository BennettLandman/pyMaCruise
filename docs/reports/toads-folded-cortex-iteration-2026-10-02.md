# Larger folded-cortex TOADS iteration checkpoint

Date: 2026-10-02  
Implementation commit: `7c593c4b6f22664320b05e3c0fcb90c0342ee1bf`  
Working directory: `pymacruise`

## What this adds

The new `tests/data/synthetic/toads_folded_cortex_iteration/` fixture is a
deterministic 27×25×23 RAS-oriented volume with a closed, folded narrow ribbon,
a one-sided notch, three tissue classes, varying memberships, and spatial
shape priors. Its continuing-iteration test composes the current implementation
through thinning, an ordinary membership update, the non-stopping control
decision, growth, and refreshed relations, boundary distances, and closest
labels.

This is useful integration coverage beyond the existing compact 13×13×11
folded slice. In the new run, the topology gate was exercised 2,094 times:
1,354 accepted and 740 rejected. The thinning queue processed 866 entries and
the growing queue processed 1,815 entries. Control selected the non-stopping
branch and enabled growing. Exact result digests and float32 control values are
stored in `python_integration_reference.json` as a regression snapshot.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_folded_cortex_iteration.py tests/test_toads_initialized_thinning_pass.py tests/test_toads_growing_pass.py tests/test_toads_object_membership_pass.py tests/test_toads_iteration_control.py tests/test_toads_relation_state_refresh.py tests/test_toads_relation_list.py tests/test_toads_relation_diffusion.py
```

Result: **23 passed in 1.46 seconds**.

The new full-volume digest is a Python integration regression reference. It is
not an output from a Java harness and should not be described as source
equivalence. Selected kernels within the composed path have their own source
comparisons: initialized thinning (`tests/test_toads_initialized_thinning_pass.py`),
growing (`tests/test_toads_growing_pass.py`), membership (`tests/test_toads_object_membership_pass.py`),
iteration decision (`tests/test_toads_iteration_control.py`), and relation
construction/refresh (`tests/test_toads_relation_list.py`,
`tests/test_toads_relation_diffusion.py`, and
`tests/test_toads_relation_state_refresh.py`). Those tests compare their
scoped behavior to selected VandyCRUISE Java harness references; they do not
validate this larger composed volume against a Java run.

## Limits and next evidence

The continuing-iteration adapter still receives fixed upstream memberships,
priors, ordering, and boundary state. This case does not cover atlas
registration, statistics dispatch, inhomogeneity fitting, stopping/class
transfer, or a complete TOADS run to convergence. It does not establish
clinical performance. The next strong checkpoint is a Java-backed intermediate
state comparison on a folded phantom, followed by connecting the actual
registered workflow to repeated TOADS outer-loop iterations.
