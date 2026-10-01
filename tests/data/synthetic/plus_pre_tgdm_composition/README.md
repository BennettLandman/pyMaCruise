# P40 pre-TGDM composition fixture

This project-generated 25×25×25 folded-cortex phantom supplies the selected
pre-TGDM composition test with asymmetric segmentation labels, T1-derived TOADS
tissue maps, multi-atlas tissue maps, and masks. All images are compressed
NIfTI (`.nii.gz`) on one anisotropic RAS+ grid (0.8, 0.8, 1.2 mm).

The case includes folded cortical interfaces and selected boundary-overlook
label 202. It is synthetic, contains no subject data, and exercises only the
bounded boundary → tissue update → ACE → GVF adapter composition. TOADS/atlas
registration and the topology-corrected level set consumed by TGDM are not
part of this fixture.

Regenerate the inputs from this directory with:

```sh
python generate_fixture.py
```

`provenance.json` records input SHA-256 digests, dimensions, and RAS affine.
