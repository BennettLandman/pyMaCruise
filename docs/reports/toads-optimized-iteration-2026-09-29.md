# P18as: one optimized-tissue TOADS iteration

## Scope

This slice composes the selected main-loop order for one optimized CSF/GM/WM
iteration: update damped centroids from the current memberships, then recompute
the image/prior terms and run the ordered membership pass with those updated
centroids. The Python API returns updated memberships, the pass objective, and
updated centroids.

The selected source call order is in `AlgorithmMultipleToads.java`: inside the
main iteration, it calls `segmentation.computeCentroids()`, then
`segmentation.computeVariances()`, then `segmentation.computeMemberships()`.
The Java harness implements those centroid-then-membership formulas on the
same deterministic synthetic 7×6×5×3 grid. The centroid calculation calls the
committed Java `PowerTable` implementation. Four JSON references cross
multiplicative field off/on with PSF ridge off/on; existing RAS `.nii.gz`
fixtures are reused.

Selected VandyCRUISE source hashes:

- `MultipleToadSegmentation.java`: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- `AlgorithmMultipleToads.java`: `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `PowerTable.java`: `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`

Regenerate the formula references with
`scripts/generate_java_toads_optimized_membership_pass_reference.sh`.

## Validation

`tests/test_toads_optimized_iteration.py` compares all four configurations.
Updated centroids match exactly. Membership maximum absolute difference is
`1.1921e-7`; objective absolute difference is at most `2.99e-8`. The tests
also check untouched outer borders. These checks validate this source-formula
slice, not full TOADS segmentation, topology admissibility, thinning,
inhomogeneity-field optimization, convergence, or atlas deformation.

The focused comparison modules pass **15 tests**. The full configured suite
passes **352 tests in 53.77 seconds** on implementation commit
`59d25c8cf8ccabaa0f3a194cb275a2645d3b0496`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md).
