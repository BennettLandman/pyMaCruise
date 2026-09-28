# LesionTOADS intensity-initialization phantom

This deterministic `9 × 8 × 7` RAS phantom contains T1 MPRAGE and FLAIR
channels with distinct modular intensity ramps, periodic negative values, and
four bright outliers. The negative values exercise the source constructor's
clamp-to-zero step. The atlas prior values are the selected lesion-atlas
entries for T1 MPRAGE and FLAIR:

- Lesion: `[0.4, 1.2]`
- Black hole: `[0.1, 0.6]`

The Java-only reference harness emits per-channel robust bounds, every
normalized voxel, and both initialized centroid vectors. Recreate the input
NIfTI files and recipe with:

```sh
.venv/bin/python tests/data/synthetic/lesion_toads_preprocessing/make_phantom.py
```

Recreate the reference artifact with:

```sh
javac -d /tmp/pymacruise-p39-reference tests/java_reference/LesionToadsPreprocessingReferenceHarness.java
java -cp /tmp/pymacruise-p39-reference LesionToadsPreprocessingReferenceHarness \
  > tests/data/synthetic/lesion_toads_preprocessing/java_reference.json
```

The fixture isolates the numeric preprocessing and centroid initialization;
it does not validate lesion atlas registration or the full topology-preserving
segmentation. The source divides by a zero robust span without checking.
PyMaCRUISE instead raises `ValueError` for that invalid input, so the
implementation intentionally differs from Java's resulting NaN/infinity
behavior in that case.
