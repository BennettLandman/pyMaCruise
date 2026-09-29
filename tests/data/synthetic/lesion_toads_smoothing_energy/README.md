# LesionTOADS six-face smoothing-energy reference

The fixture contains four selected-Java one-voxel cases for `OBJECT`,
`LESION`, `BLACKHOLE`, and `OUTMASK` classes. Six synthetic face-neighbor
membership vectors exercise source face traversal order, current-class
exclusion from the other-class sum, lesion and black-hole auxiliary
cross-coupling, the lesion-type self-neighbor term, and full versus half
energy updates.

No iterative smoothing or membership change is performed. Regenerate with
`scripts/generate_java_lesion_toads_smoothing_energy_reference.sh`.
