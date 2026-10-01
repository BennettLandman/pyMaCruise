# Initial-membership small-distance threshold fixture

This is a deterministic 4×4×4 two-channel RAS NIfTI with a 2×2×2 active crop
and one zero padding voxel on each face. The values intentionally make class
squared distances exactly zero, nonzero but below `LesionToadPreprocess.ZERO`
(`1e-20f`), and above it. It isolates the source's reciprocal-versus-`INF`
decision within the existing fixed-input initial membership method.

`generate_fixture.py` creates the NIfTI and centroid metadata. The Java
reference harness is pinned to VandyCRUISE
`LesionToadPreprocess.java` SHA-256
`72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08` and
compares complete output membership bits. It does not extend scope to crop
creation, preprocessing, centroid estimation, or atlas registration.
