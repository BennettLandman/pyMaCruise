# Bundled TOADS atlas RAS loading smoke check

Code tested: `845eecfa43b9862ce81e7d9322aeedccf4593b7b`.

The selected bundled descriptor
`data/atlases/cruise-atlas-10obj-toads2008.txt` loaded through
`load_toads_atlas`, `load_toads_atlas_ras_data`,
`load_toads_intensity_priors`, and
`toads_atlas_shapes_in_structure_order`. Assertions verified topology grid
`(172, 216, 164)`, 11 descriptor-ordered shape classes, finite T1_MPRAGE prior
shape `(1, 11)`, and topology labels `{1, 5, 6, 14, 15, 16, 17, 18, 22, 24, 25}`.
The template spacing returned from the descriptor is `(0.9375, 0.9375, 0.9375)` mm.
The loaded topology/maps and stacked float32 class maps required 542,270,592
array bytes.

The smoke command completed successfully and printed:

```text
atlas=cruise-atlas-10obj-toads2008.txt dimensions=(172, 216, 164) classes=11 priors=(1, 11) topology_labels=[1, 5, 6, 14, 15, 16, 17, 18, 22, 24, 25]
voxel_spacing_mm=(0.9375, 0.9375, 0.9375) array_bytes=542270592
```

This checks bundled asset decoding, descriptor order, labels, and prior parsing.
It is not a pytest run, reads no subject MRI, and provides no segmentation or
registration accuracy evidence.
