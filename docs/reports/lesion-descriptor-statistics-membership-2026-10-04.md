# LesionTOADS descriptor statistics-to-membership transition

`run_lesion_toads_descriptor_statistics_membership` composes the selected
post-thinning stages for the bundled 13-class descriptor: refresh class and
auxiliary centroids, variances and channel weights; calculate centroid
similarities; compute active raw class terms; run the descriptor-aware
membership pass; and calculate the source energy ratio.

The selected `use_lesion_weight` option routes the descriptor's FLAIR channel
and per-channel lesion PowerTable exponent through both statistics refresh
and raw lesion energy. Missing exponents and out-of-range channels are
rejected. Tests compare the orchestration with the validated Python
components, exercise weighted and unweighted modes, check finite energy
ratios and verify that input arrays are not mutated.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_active_class_energy.py tests/test_lesion_descriptor_statistics_membership.py tests/test_lesion_data_energy_maps_special.py
```

Focused result: **8 passed, 0 failed, 0 skipped in 0.69s**.

Full suite command: `./.venv/bin/pytest -q`.

Full suite result: **1,153 passed, 0 failed, 0 skipped, 11 expected RAS
warnings in 828.21s** on tested commit
`936ba20da8ea6ab3ea3bdb238cc1ac100790c0dc`. Peak RSS was not sampled. See the
[full-suite report](full-pytest-2026-10-04-descriptor-statistics-membership.md).

This is synthetic, fixed-state integration and uses the real bundled atlas
descriptor for class order. It does not provide a new Java harness comparison
for the complete adapter. Thinning, topology/relation rebuilding, correction
field updates, lesion-distance modulation, repeated outer iterations,
convergence, complete LesionTOADS execution, MaCRUISEplus, and real MRI remain
outside this row.
