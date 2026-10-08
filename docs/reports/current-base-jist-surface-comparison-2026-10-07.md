# JIST surface comparison on the brain-resolution synthetic base core

This compares selected Java `IsoSurfaceOnGrid.solveOriginal(..., 18, 0, false)` with the production Python Lewiner extractor on the exact same committed Java TGDM outputs. The inputs are synthetic, RAS, and `(172, 216, 164)` at 0.9375 mm spacing; this is a stage comparison, not an end-to-end Java run.

JIST source SHA-256: `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`. Wrapper SHA-256: `3f7b54e099f2550f42a572e96070642b8124a93ab1b45eca1aa544ecec1758b4`.

| Surface | JIST V/F/χ | Lewiner V/F/χ | Symmetric sampled vertex-to-triangle mean (mm) | Sampled max (mm) | Matched unordered faces |
|---|---:|---:|---:|---:|---:|
| inner | 116166/232328/2 | 116178/232368/-6 | 0.00000837 | 0.438229 | 67.94% |
| central | 121074/242144/2 | 121100/242224/-12 | 0.00001699 | 0.403988 | 63.79% |
| outer | 123828/247652/2 | 123872/247828/-42 | 0.00003110 | 0.425185 | 67.52% |

The symmetric sampled means are below 0.01 mm and sampled maxima are below 1 mm for all three surfaces, meeting the owner-provided P08 thresholds on these fixtures. The Lewiner topology counts differ: its mesh has multiple connected components and nonmatching Euler characteristics where JIST produces one connected genus-zero surface. Only about 63.8–67.9% of JIST triangles map after nearest-vertex matching. These remain material connectivity differences under P37 even though the sampled geometric position metrics pass.

Applying the same Python closest-point thickness operator to each nested mesh set gives JIST mean `1.649141538 mm` and Lewiner mean `1.649053461 mm`, a delta of `-0.000088077 mm`. P08 has no thickness-specific tolerance. This is a geometry-impact measurement, not Java thickness equivalence.

The Java harness compiles the selected source body directly, with minimal image/progress and normal-generation stubs. It writes mesh arrays under ignored `build/`; this report stores only hashes and metrics. Vertex-to-triangle distances sample vertices and are not a continuous Hausdorff certificate. Reproduce with `./.venv/bin/python scripts/compare_jist_surface_current_base.py`; raw results: [JSON report](current-base-jist-surface-comparison-2026-10-07.json).

Synthetic selected-source evidence only. It does not establish complete MaCRUISE workflow equivalence or validate human MRI outputs.
