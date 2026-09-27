# Segmentation refinement synthetic fixture

Three committed `.nii.gz` files provide a compact 11³ RAS grid with 0.5, 0.75,
and 1.0 mm voxel spacings. The integer MaSeg labels and inner/outer level sets
place known targets for inner WM filling, modal-label ties, outer label removal,
central WM/GM filling, an unresolved border case, and preservation of an
unrelated label. No subject data are present.

Regenerate with `python tests/data/synthetic/refinement/make_phantom.py` from
the repository root. Verify with `shasum -a 256 -c SHA256SUMS`.
