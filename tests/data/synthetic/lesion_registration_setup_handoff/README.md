# LesionTOADS final initialization registration handoff

This is an all-zero 5 x 6 x 7 x 3 float32 membership buffer, matching the
segmentation constructor state before the wrapper computes initial
memberships. It is stored in compressed NIfTI on an RAS grid with 0.8 x 1.0 x
1.2 mm voxels.

The selected default wrapper passes this membership reference to
`LesionToadAtlas.updateShapeRegistration` with rigid mode, one iteration, one
level, and coupling 1.0; that method fixes subsample to 2. This fixture checks
only the wrapper handoff fields and array reference. It does not claim to port
the atlas transform-model reconfiguration or optimize a registration.

The fixture includes selected source excerpts for wrapper defaults and call
order, constructor allocation/getter behavior for the membership buffer, and
the initial field assignments in `LesionToadAtlas.updateShapeRegistration`.
The test pins each excerpt hash and checks that the source fields and call
sequence agree with the harness values. The NIfTI fixture SHA-256 is also
pinned so the reference remains tied to the committed zero membership map.
