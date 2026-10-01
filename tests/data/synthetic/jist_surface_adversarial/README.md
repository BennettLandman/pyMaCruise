# Adversarial JIST-versus-Lewiner synthetic surfaces

The comparison runner generates four float32 level sets on a 41×41×41,
1 mm identity-affine RAS grid:

- a genus-one torus;
- a two-lobe dumbbell joined by a 3 mm connector;
- two radius-7 mm components with a 0.2 mm continuous gap;
- a sphere centered at x=2 mm whose surface crosses the image boundary.

The selected Java algorithm is `IsoSurfaceOnGrid.solveOriginal` with
connectivity 18, iso-level 0, and `inclusive=false`, the settings traced in
`GetMeshSurfaceFromLevelSet`. The runner writes the synthetic `.nii.gz`
inputs and compares both implementations on those same arrays. It records
connected components, Euler characteristic, boundary/nonmanifold edges,
nearest-vertex mapped triangle agreement, and symmetric vertex-to-triangle
sampled distances. Java failures and Python failures are recorded per case so
one boundary failure does not hide other comparisons.

Run from `pymacruise` with:

```sh
.venv/bin/python scripts/compare_jist_adversarial_surfaces.py
```

The report is `docs/reports/surface-extraction-adversarial-java-comparison-2026-09-30.md`;
the machine-readable results and source/input checksums are in
`comparison.json`. Raw arrays, extracted Java meshes, and class files are
written under ignored `build/` only.

These measurements are substitution-risk evidence. Vertex-sampled
vertex-to-triangle distances are not a continuous Hausdorff bound, and a
matching metric on these phantoms is not P37 acceptance or end-to-end
equivalence.
