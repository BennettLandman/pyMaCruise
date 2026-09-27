# Java TGDM golden images

These generated, synthetic-only reference volumes are produced by
`scripts/generate_java_tgdm_reference.sh`. They use a 25×25×25 phantom, 1 mm
isotropic voxels, identity RAS affine, and float32 NIfTI data. The Java source
and fixture-generation details are documented in
[`tests/java_reference/README.md`](../../../../java_reference/README.md).

Each `.nii.gz` file records one source TGDM output stage. No subject data is
included. `SHA256SUMS` covers the compressed images in this directory.
