# Terminal synthetic PLUS surface extraction: Java JIST and Python source port

The source-selected `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)` and the production Python source-derived 18/6 case-template extractor were given the exact same terminal checkpoint level-set arrays. This is a synthetic stage comparison of the current source port; it does not establish upstream or whole-workflow equivalence.

The registered LesionTOADS loop stopped by `convergence` at iteration 10. The input level sets were checked by SHA-256 before comparing the generated meshes.

| Surface | Java JIST V/F/components/χ | Python source-port V/F/components/χ | Symmetric sampled mean (mm) | Sampled max (mm) |
|---|---:|---:|---:|---:|
| inner | 128084/256184/4/-8 | 128084/256184/4/-8 | 0.00000001 | 0.000244 |
| central | 142532/285080/4/-8 | 142532/285080/4/-8 | 0.00000003 | 0.000244 |
| outer | 197562/395140/4/-8 | 197562/395140/4/-8 | 0.00000001 | 0.000173 |

Using the same Python closest-point thickness operator on both mesh sets, mean thickness is 4.141158136 mm (JIST) and 4.141157975 mm (Python), delta -0.000000161 mm. This is not a Java thickness comparison.

Distance summaries sample mesh vertices against triangles; their maxima are not continuous Hausdorff certificates. See the adjacent JSON for hashes, topology, and method provenance.
