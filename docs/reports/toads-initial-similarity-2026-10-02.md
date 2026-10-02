# TOADS initial centroid similarity (P18cyap)

P18cyap carries the selected base-TOADS startup through the centroid-similarity
prelude inside `computeExactMemberships()`. It takes the initialized centroids,
assigned-class standard deviations, and image-wide inverse-variance weights,
then produces the descriptor-order pairwise matrix with the selected default
`priorScale=0.05`. The matrix includes MASK and optimized tissue classes; its
diagonal is zero as in the Java loop. The adapter explicitly limits this stage
to the base non-lesion descriptor types. Lesion-specific pair behavior and the
voxelwise exact membership/energy sweep remain outside this stage.

Focused synthetic validation:

- Command: `./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_optimized_similarity.py tests/test_toads_pair_similarity.py`
- Result: **8 passed, 0 failed, 0 skipped in 7.82 seconds**.
- Implementation commit: `8663b867a912c5a5ff7a8ac0b009b75512e0fd15`.

The startup test verifies call order, matrix shape, diagonal, and finite values
on a compact two-class RAS atlas case. The optimized-similarity test compares
the numerical kernel to its pinned Java-reference matrix and checks source
IEEE behavior for zero deviations. This validates a two-class synthetic
composition plus the source-checked kernel; it is not an 11-class startup
comparison or a complete TOADS membership update.
