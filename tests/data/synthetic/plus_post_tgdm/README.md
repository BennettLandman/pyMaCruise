# MaCRUISEplus post-TGDM synthetic inputs

These six 13³ NIfTI images are fixed cropped outputs representing the inputs
to the post-TGDM part of the selected MaCRUISEplus pipeline. They are derived
from the committed 17³ anisotropic RAS nested-sphere level sets in
`../plus_surfaces/`, cropped to voxel bounds `[2, 14]` on each axis. The
fixture uses 1.0 × 1.5 × 2.0 mm spacing and a nonzero crop origin; all images
are RAS-oriented.

The integration test starts from these files and composes crop restoration,
nested surface extraction, and thickness measurement. It does not execute or
claim TOADS/LesionTOADS, MaACE/GVF, or full MaCRUISEplus. Surface extraction
uses scikit-image Lewiner marching cubes; equivalence with JIST
`IsoSurfaceOnGrid` remains an explicit unresolved risk.

Regenerate the files by running `make_phantom.py` from the repository root's
`pymacruise` directory. `SHA256SUMS` records their committed content.
