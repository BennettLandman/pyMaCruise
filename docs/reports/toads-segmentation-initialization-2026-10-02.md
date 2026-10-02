# TOADS segmentation initialization state (P18cyak)

**Implementation/test commit:** `04e7e9b346cbca0ae7285e4cf2735f6e3c56d7b7`  
**Validation date:** 2026-10-02

This slice composes the selected `AlgorithmMultipleToads` initialization after
the second atlas registration. It copies regular centroids into both the
segmentation centroid and prior-centroid state, carries the preprocessing
robust maximum and scale, maps the descriptor lesion prior through the
source-normalized intensity formula, and computes the fresh exact membership
array. The selected base wrapper sets `useApprox=false` and `approx=classes`,
so the source atlas-aware overload dispatches to exact membership
initialization. `load_toads_lesion_priors` follows the descriptor `Lesions`
section and returns zeros when it is absent or a requested modality has no
listed value.

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_atlas_priors.py tests/test_toads_initialization.py tests/test_toads_lesion_centroid_initialization.py tests/test_toads_alignment_startup.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_registration_sampling.py tests/test_toads_registration.py
```

Result: **59 passed, 0 failed, 0 skipped in 6.93 seconds** on the commit above.
The compact synthetic atlas verifies call order, source defaults, centroid and
intensity state, zero-prior fallback, lesion-prior aliases, non-finite
lesion-prior rejection, and equality of the fresh membership array with the exact memberships previously
used by the second registration. The centroids, lesion-prior mapping, and
membership kernels have separate selected-source comparisons. This combined
state handoff is synthetic composition evidence; it does not establish a
complete TOADS iteration or Java end-to-end segmentation equivalence.
