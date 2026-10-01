# P18czo synthetic default-atlas approximate-membership reference

This fixture validates the atlas-aware 11-class `initialApproxMemberships`
method only. `images.nii.gz` and `transformed_shapes.nii.gz` are asymmetric
7×6×5 RAS volumes with two image channels and the selected TOADS2008 class
order. The atlas descriptor determines the class order. The first rigid
registration pyramid calls the separate no-atlas overload (P18czp); the later
atlas-aware call is the one composed with registration in P18czt.

The Java harness pins every initial sparse membership and candidate slot for
the full fixture, plus shape-seeded ties, zero-cost reciprocal fallback, and a
small center-only process-mask adapter case. The selected Java method itself
does not accept a boolean process mask; the adapter case verifies that
unprocessed Python voxels stay zero. Input identity hashes for the harness,
Numerics source, NIfTI files, and atlas descriptor are recorded in
`java_reference.json` and checked by the focused test.

Regenerate the reference with:

```sh
sh tests/java_reference/generate_toads_default_approximate_membership_reference.sh
```

This reference does not validate registration, later exact-membership sweeps,
or complete TOADS segmentation.
