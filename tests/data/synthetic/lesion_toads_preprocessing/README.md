# LesionTOADS intensity-initialization phantom

This deterministic `9 × 8 × 7` RAS phantom contains T1 MPRAGE and FLAIR
channels with distinct modular intensity ramps, periodic negative values, and
four bright outliers. The negative values exercise the source constructor's
clamp-to-zero step. The atlas prior values are the selected lesion-atlas
entries for T1 MPRAGE and FLAIR:

- Lesion: `[0.4, 1.2]`
- Black hole: `[0.1, 0.6]`

The selected `AlgorithmLesionToads` source pins `robustMaxMin=false`. Its call
order is constructor negative clamp and ordinary min/max, `normalizeImages`
using `(image - Imin) / (Imax - Imin)`, crop, tissue-centroid initialization,
then the separate lesion and black-hole centroid methods. Those latter
methods compute robust bounds on the normalized (in the wrapper, cropped)
channels and use `low + prior * (high - low)`. The fixture harness reproduces
normalization and those centroid formulas on the full grid; it does not
reproduce cropping. It emits ordinary bounds, centroid robust bounds, every
normalized voxel, and both initialized centroid vectors. Source SHA-256 values
and the wrapper call order are pinned in `java_reference.json`.

Recreate the input
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

The fixture isolates the active normalization and centroid formulas; it does
not validate cropping, lesion atlas registration, or the full
topology-preserving segmentation. The active source divides by a zero
ordinary or robust centroid span without checking. PyMaCRUISE instead raises
`ValueError` for those invalid inputs, so it intentionally differs from
Java's resulting NaN/infinity behavior. The robust normalization and
`low / span + prior` formulas in the Java class are an inactive
`robustMaxMin=true` branch and are not represented as selected-wrapper behavior.
