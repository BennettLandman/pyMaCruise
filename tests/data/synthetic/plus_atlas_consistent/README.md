# Atlas-aligned synthetic PLUS input case

This fixture is generated on the bundled LesionTOADS topology atlas grid. The
atlas topology IDs use the source-derived MaSeg grouping: topology 15 alone is
cortical GM (100), 25–27 are cortical WM (44), 5–6 are CSF (5), topology 1 is
background (0), and the remaining labeled structures are other tissue (75).
T1 MPRAGE and FLAIR intensities use the matching
LesionTOADS descriptor means with a small deterministic smooth within-class
modulation. FLAIR uses an exterior intensity baseline of 80 to avoid negative
values from cubic registration interpolation at the brain boundary. The brain
mask is the non-background atlas support. All inputs are RAS `.nii.gz` files
on a shared grid.

This is intended to probe the workflow when labels, intensities, and atlas
priors have aligned anatomy. Because it is derived from the same atlas, it is
not independent ground truth for atlas registration, TOADS, or whole-workflow
source equivalence.

Regenerate with:

```sh
PYTHONPATH=src .venv/bin/python tests/data/synthetic/plus_atlas_consistent/generate_fixture.py
```
