# TOADS constructor and relation-state composition (P18cyal)

**Implementation commit:** `cd8d9bdf1854e2c0e79118780071c0868de7c8e1`  
**Validation date:** 2026-10-02

The atlas-backed synthetic startup now carries the refreshed RAS topology into
the selected segmentation constructor state. Descriptor topology labels map
to zero-based class indices in descriptor order; the adapter initializes the
segmentation, skeleton and relation maps, zero availability and boundary
distances, and per-class volumes. It then calls the already source-checked
relation, boundary-distance and four-slot closest-label adapters in source
order with the selected limit of 3.

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_initial_relation_state.py tests/test_toads_relation_list.py tests/test_toads_boundary_distances.py tests/test_toads_closest_labels.py tests/test_toads_relation_state_refresh.py
```

Result: **32 passed, 0 failed, 0 skipped in 9.42 seconds** on the implementation
commit above. Tests verify descriptor-order mapping, seed arrays, volumes,
unknown-label rejection, and relation-state dimensions and outputs on the
compact asymmetric synthetic atlas. The component kernels have selected Java
reference coverage; this composition is synthetic wiring evidence and is not
a direct Java comparison of the complete segmentation constructor. Variance,
inhomogeneity-field setup and the initial membership/energy pass remain open.
