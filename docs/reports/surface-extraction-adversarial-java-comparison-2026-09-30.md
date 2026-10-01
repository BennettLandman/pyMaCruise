# JIST `IsoSurfaceOnGrid` versus Lewiner: adversarial synthetic level sets

Selected source settings: connectivity 18, iso-level 0, `inclusive=false`; all inputs are identical float32 arrays on a 1 mm identity RAS affine.

| Shape (expected topology) | JIST V/F/χ; components; boundary / nonmanifold edges | Lewiner V/F/χ; components; boundary / nonmanifold edges | Mapped faces | Symmetric sampled vertex→triangle mean / max (mm) | Result |
|---|---|---|---:|---:|---|
| torus (1 comp, χ=0, closed) | 2312/4624/0; 1; 0 / 0 | 2312/4624/0; 1; 0 / 0 | 3136/4624 (67.8%) | 1.88148311e-07 / 1.63106643e-06 | both ran |
| narrow_neck_dumbbell (1 comp, χ=2, closed) | 2486/4968/2; 1; 0 / 0 | 2486/4968/2; 1; 0 / 0 | 3568/4968 (71.8%) | 1.70637764e-07 / 1.30304056e-06 | both ran |
| near_touching_components (2 comp, χ=4, closed) | 1764/3520/4; 2; 0 / 0 | 1764/3520/4; 2; 0 / 0 | 2512/3520 (71.4%) | 1.75603029e-07 / 1.52587891e-06 | both ran |
| boundary_crossing_sphere (1 comp, χ=1, open clipped patch) | 945/1820/1; 1; 68 / 0 | 945/1820/1; 1; 68 / 0 | 1248/1820 (68.6%) | 1.61540924e-07 / 7.64123191e-07 | both ran |

Both methods matched the expected components and Euler characteristic for all four phantoms, with zero nonmanifold edges. The three interior surfaces are closed; the boundary-crossing sphere produced the expected open clipped patch with 68 boundary edges in each method. A sampled distance or topology match does not certify continuous surface Hausdorff distance. A difference in mapped faces is evidence of triangulation/connectivity sensitivity in the Lewiner substitution.

Source/input metrics and per-case failures are in [`comparison.json`](../../tests/data/synthetic/jist_surface_adversarial/comparison.json); generator and runner: [`compare_jist_adversarial_surfaces.py`](../../scripts/compare_jist_adversarial_surfaces.py). Selected source SHA-256: `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`; wrapper SHA-256: `3f7b54e099f2550f42a572e96070642b8124a93ab1b45eca1aa544ecec1758b4`.

The reference runner writes raw inputs, generated Java meshes, and class files only under ignored `build/`; committed fixture inputs are the four synthetic `.nii.gz` volumes.
