# Mixed LesionTOADS centroid phantom

`mixed_centroid_inputs.nii.gz` is a committed synthetic RAS NIfTI with a
6 × 5 × 4 grid and 21 frames: hard labels, 13 class memberships, two image
channels, two bias fields, an OUTLIER map, LESION map, and BLACKHOLE map. It
contains all selected `computeCentroids()` categories, labels every class
except class 12, links categories across centroid groups, and exercises
auxiliary winner ties/routes plus a below-ZERO optimized factor.

The fixture is created by
`scripts/generate_lesion_mixed_centroid_fixture.py`. The selected Java
`computeCentroids()` behavior and the exact upstream MIPAV `PowerTable` class
are exercised by `LesionMixedCentroidsReferenceHarness.java`; regenerate the
committed raw-float-bit reference with
`scripts/generate_java_lesion_mixed_centroids_reference.sh`. The source file
SHA-256 is `77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0`;
the upstream `PowerTable.java` SHA-256 is
`2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`.

These phantom checks establish agreement with the selected Java harness, not
equivalence of the full historical pipeline or validation on acquired scans.
