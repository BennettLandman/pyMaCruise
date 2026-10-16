# AtlasEM Gaussian membership update slice — 2026-10-16

## Source-backed behavior

`update_atlas_em_gaussian_memberships` ports the selected one-channel, no-field, Atlas-prior branch of `AtlasEMSegmentation.computeMemberships` from source class SHA-256 `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a`.

The fixed-input slice includes scalar covariance inversion, Gaussian log likelihood, atlas-prior weighting, posterior normalization with the source `1e-30f` denominator seed, source-order x/y/z/class iteration, in-place six-neighbor exclusion smoothing, clearing masked voxels, cumulative energy, and maximum membership change. It requires the false six-face padding that the source loop relies on for neighbor indexing. The Atlas prior and covariance remain caller-supplied; this code does not register or evolve the atlas.

`tests/java_reference/AtlasEMGaussianMembershipReferenceHarness.java` copies this selected 1-channel branch into a dependency-free Java harness. On a deterministic 4×4×4 three-class phantom with positive priors in the 2×2×2 interior and smoothing `0.625`, all output memberships and both metrics match float32 bits exactly. The source excerpt does not exercise zero-prior / zero-old-membership NaN energy cases; the Python function documents and retains Java's `log(0) * 0` NaN behavior rather than suppressing it.

Focused validation on the working tree based on `5abcef6b`:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_outputs.py tests/test_mipav_morphology.py
```

Result: **12 passed in 2.87 s**.

## Limits

This is one fixed-input EM update and not AtlasEM registration, deformable atlas evolution, multi-channel/inhomogeneity-field membership updates, stopping logic, classification output, the pre-LesionTOADS Dura stage, or full PLUS integration. The selected `ibsr8obj` atlas descriptor and payload are unavailable, so this slice is not evidence of a runnable source-configured AtlasEM stage.
