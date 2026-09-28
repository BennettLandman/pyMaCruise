# Registration transform phantom

This committed synthetic fixture tests transform direction, label and
continuous interpolation, point/surface coordinates, fixed-space masks, and
inverse transforms. The moving landmark is RAS+ `(5, 8, 7) mm`; point
application maps it to fixed RAS+ `(7, 8, 7) mm`. ANTs image resampling uses
the opposite fixed-grid-to-moving-image sampling direction, so its forward
image transform stores an LPS translation of `(2, 0, 0) mm`. Images are 17³,
1 mm, and RAS.

Regenerate with `python tests/data/synthetic/registration/make_phantom.py`
from the repository root. Verify fixture hashes using
`shasum -a 256 -c SHA256SUMS`.
