# Brain-resolution synthetic FLAIR registration phantom

This committed `(172, 216, 164)` identity-affine RAS fixture is analytic
synthetic data, not a human scan. The T1 image contains nested ellipsoidal
tissue regions with a folded cortical modulation. The FLAIR image uses a
monotonic contrast remapping, one synthetic focal hyperintensity, and a known
integer-voxel shift relative to T1. It exercises the rigid checkpoint at
brain-resolution dimensions. It does not model clinical noise, anatomy, or
prove NiftyReg/Aladin equivalence.

Regenerate both `.nii.gz` files with:

```bash
./.venv/bin/python tests/data/synthetic/plus_flair_brain_resolution/generate_fixture.py
```
