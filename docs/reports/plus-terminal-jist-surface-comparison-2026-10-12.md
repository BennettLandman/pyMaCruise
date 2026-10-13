# JIST surface extraction on the corrected terminal synthetic PLUS case

The source-selected `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)` and Python Lewiner extractor were given the exact same terminal checkpoint level-set arrays. This is a synthetic stage comparison; it isolates the P37 triangulation substitution and does not establish whole-workflow equivalence.

The registered LesionTOADS loop stopped by `convergence` at iteration 10. The input level sets were checked by SHA-256 before comparing the generated meshes.

| Surface | JIST V/F/components/χ | Lewiner V/F/components/χ | Symmetric sampled mean (mm) | Sampled max (mm) |
|---|---:|---:|---:|---:|
| inner | 128084/256184/4/-8 | 128084/256184/4/-8 | 0.00000001 | 0.000244 |
| central | 142532/285080/4/-8 | 142532/285080/4/-8 | 0.00000003 | 0.000244 |
| outer | 197562/395140/4/-8 | 197562/395140/4/-8 | 0.00000001 | 0.000173 |

Using the same Python closest-point thickness operator on both mesh sets, mean thickness is 4.141158 mm (JIST) and 4.141158 mm (Lewiner), delta -0.000000 mm. These are geometry-isolation measurements, not Java thickness equivalence.

Distance summaries sample mesh vertices against triangles; their maxima are not continuous Hausdorff certificates. See the adjacent JSON for hashes, topology, and method provenance.
