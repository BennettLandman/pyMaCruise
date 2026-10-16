# AtlasEM preprocessing source slices — 2026-10-16

## Evidence

`src/pymacruise/atlas_em_preprocessing.py` implements selected initialization behavior from `AtlasEMPreprocess.java`, whose SHA-256 is `acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b`.

- Positive-only image clipping and per-channel robust min/max are delegated to `robust_intensity_bounds`, already source-compared for `ImageFunctions.robustMinimum/robustMaximum` by the TOADS intensity harness.
- Normalization uses the source robust span `Ihigh-Ilow` and preserves zero background; it does not subtract `Ilow`.
- Active `robustMaxMin=true` centroid initialization computes `Ilow/(Ihigh-Ilow)+prior` in float32.
- Initial memberships reproduce inverse squared Euclidean distance, `ZERO=1e-20`, `INF=1e20`, exact first-channel centroid grouping, and source float32 operation order. The default exclusive bounds leave the last index of each spatial axis zero, matching `xNi=nix-1` and analogous Y/Z source loops; callers can provide explicit cropped interior bounds.

`tests/java_reference/AtlasEMPreprocessReferenceHarness.java` contains excerpts of the selected centroid and membership formulas. `tests/test_atlas_em_preprocessing.py` compares every initialized centroid and every membership float bit-for-bit against that harness on a deterministic 4×3×3 two-channel phantom. It checks default crop bounds and tied centroids. A separate test covers robust-span normalization and degenerate-span rejection.

Focused validation command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_preprocessing.py tests/test_atlas_em_statistics.py \
  tests/test_atlas_em_outputs.py tests/test_mipav_morphology.py
```

Result on the working tree based on `2ed21b2e`: **10 passed in 2.54 s**.

## Limits

This is not the complete AtlasEM preprocessing object or workflow. Image cropping, padding/index mapping, atlas resampling and registration, topology initialization, AtlasEM iterations and stopping criteria, and the pre-LesionTOADS Dura evolution are still open. The selected `ibsr8obj` descriptor and referenced atlas payload are missing, so the initialization cannot yet be assembled into a source-configured end-to-end PLUS run. The Java harness is a pinned method excerpt, not an independent MIPAV/JIST execution.
