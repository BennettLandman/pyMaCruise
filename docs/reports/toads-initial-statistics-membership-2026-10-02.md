# TOADS initial statistics and exact membership composition (P18cyaq)

P18cyaq adds a source-ordered composition over an explicit initialized
11-class state. It computes assigned-class deviations, full-grid
image-variance weights, and the default-scale pairwise centroid-similarity
matrix, then passes those results to one exact membership/energy sweep.
Transformed shapes, relation maps, closest labels, intensity exponents,
modality weights, and atlas/smoothing priors remain fixed upstream inputs.

Validation uses the asymmetric folded 9×8×7 two-channel RAS fixture from the
selected 11-class dispatcher. All eleven assigned classes are represented in
the synthetic segmentation label state. The test checks stage order, output
dimensions/finiteness, and that the composed pass matches a direct call with
the computed statistics and similarity. It does not claim a combined Java
wrapper oracle; the individual numerical kernels and exact dispatcher have
separate pinned selected-Java reference tests.

- Command: `./.venv/bin/pytest -q tests/test_toads_initial_statistics_membership.py tests/test_toads_default_atlas_dispatch.py tests/test_toads_optimized_class_variance.py tests/test_toads_optimized_variance.py tests/test_toads_optimized_similarity.py`
- Result: **18 passed, 0 failed, 0 skipped in 4.40 seconds**.
- Implementation commit: `6fd8b12b6b048f3cf8f25bcb041b0d073542e13b`.

Atlas transform sampling, field fitting, and iterative TOADS evolution remain
outside this composition.
