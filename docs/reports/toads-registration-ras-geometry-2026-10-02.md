# P18cyag: RAS geometry in the TOADS registration setup

Implementation commit: `f6b8848bea1a3d76e9ea1d2c3e1406bc28b5bce8`.

The setup record now exposes the inputs corresponding to the selected Java
`setImageInfo` and rigid `ParametricTransform`: subject membership-grid size,
template-prior grid size, MIPAV image orientation (`AXIAL=0`), MIPAV RAS axis
enums `(L2R=2, P2A=3, I2S=5)`, orientation signs `(-1,-1,+1)`, signed subject
resolutions, descriptor template resolutions, and centers at half the voxel
counts. The setup geometry is built from the startup and the actual arrays
passed to registration, avoiding inferred grid metadata.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_registration_sampling.py tests/test_toads_registration.py
```

Result on the exact implementation commit: **43 passed, 0 failed, 0 skipped in
1.22 seconds**. The folded synthetic adapter test checks dimensions, orientation
enums/signs, signed spacing, template spacing, and both centers. The same group
includes selected Java-reference signed landmark and registration-kernel checks.

No actual subject MRI was used. This records and checks geometry handoff; it
does not establish full wrapper-state or workflow equivalence.
