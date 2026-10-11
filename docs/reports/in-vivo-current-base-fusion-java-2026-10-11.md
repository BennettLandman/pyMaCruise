# Corrected human base tissue fusion: held-fixed Java comparison

## Result

The selected Java `CombineMAandTOADS.execute()` path was run on the six
cropped inputs saved by the corrected human base preparation stage. The
fixed-input Java and Python results match bit-for-bit for GM, WM fill, and
mask across all **3,195,483 voxels per output**. Each output has zero changed
voxels and zero absolute error. Inputs and outputs are on the same RAS grid
with shape `(141, 173, 131)`.

The corrected preparation stage includes MaSeg CSF-boundary construction and
cleanup around fusion. Those extra stages are accepted by the comparator, but
this result validates only the fusion call: MaSeg cortical separation,
nearest-neighbor transfer, shared cropping, and whole-workflow equivalence
are not established here. The Java harness uses the selected 3D center plus
six-neighbor erosion and GM confident level `0.001`.

## Reproduction

```sh
.venv/bin/python scripts/compare_base_tissue_preparation_saved_java.py \
  --preparation-dir build/invivo-base-source-v3-human-boundary-2026-10-11/stages/base-tissue-preparation/attempt-01 \
  --output-dir build/invivo-base-source-v3-human-boundary-2026-10-11/java-fusion-held-fixed
```

The Java harness is capped at `-Xmx768m` and a 2 GiB sampled-RSS watchdog.
The recorded Java peak was **0.242 GiB**. The machine-readable comparison
records the exact input and provenance hashes at
`build/invivo-base-source-v3-human-boundary-2026-10-11/java-fusion-held-fixed/comparison.json`.
The human images and comparison intermediates remain in ignored `build/`.

This is source-equivalence evidence for one held-fixed fusion stage, not
evidence that the human reconstruction passes anatomical acceptance.
